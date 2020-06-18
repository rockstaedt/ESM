using JuMP
using Clp
using Plots, StatsPlots
using DataFrames, CSV


# Preprocessing

### data load ###
datapath = "data"
tspath = joinpath(datapath, "time_series")

technologies = joinpath(datapath,"technologies.csv") |> CSV.read
zones = joinpath(datapath,"zones.csv") |> CSV.read
lines = joinpath(datapath,"ntc.csv") |> CSV.read

timeseries = Dict(splitext(files)[1] => CSV.read(joinpath(tspath, files))
    for files in readdir(tspath))

######################################################################
### create sets based on input data ###

### sets ###
# |> Funktion wird angewendet
T = 1:size(timeseries["load"], 1) |> collect
TECH = technologies[:,:technology] |> Vector
DISP = technologies[technologies[:,:dispatchable] .== 1, :technology]
NONDISP = technologies[technologies[:,:dispatchable] .== 0, :technology]
S = technologies[technologies[:,:storage_efficiency] .!= 0, :technology]
CHP = technologies[technologies[:,:heat_ratio] .!= 0, :technology]
Z = unique(zones[:,:zone])
P = zones[:,:id] |> Vector


### parameters ###
### utility functions and mappings ###
zipcols(df::DataFrame, x::Symbol) = df[:,x] |> Vector
zipcols(df::DataFrame, x::Vector) = zip(eachcol(df[:,x])...) |> collect

function dictzip(df::DataFrame, x::Pair)
    dictkeys = zipcols(df, x[1])
    dictvalues = zipcols(df, x[2])
    return zip(dictkeys, dictvalues) |> collect |> Dict
end

map_id2tech = dictzip(zones, :id => :technology)
map_id2z = dictzip(zones, :id => :zone)

P_DISP = filter(x-> map_id2tech[x] in DISP, P)
P_NONDISP = filter(x-> map_id2tech[x] in NONDISP, P)
P_S = filter(x-> map_id2tech[x] in S, P)
P_CHP = filter(x-> map_id2tech[x] in CHP, P)

map_zone2plist = Dict(z => filter(x-> map_id2z[x] == z, P) for z in Z)

### core parameters ###
mc = dictzip(technologies, :technology => :marginal_cost)
g_max = dictzip(zones, :id => :g_max)

eff = dictzip(technologies, :technology => :storage_efficiency)
stor_max = dictzip(zones, :id => :storage_capacity)

heat_ratio = dictzip(technologies, :technology => :heat_ratio)

### time series ###
coldict(df::DataFrame) = Dict(string(name) => Vector(vec) for (name,vec) in pairs(eachcol(df)))

elec_demand = coldict(timeseries["load"])
heat_demand = coldict(timeseries["heat"])

g_max_by_z_tech = Dict((map_id2z[id],map_id2tech[id]) => g_max[id] for id in P)

feed_in_by_z_nondisp = Dict(
    (z,ndisp) => timeseries[ndisp][!,z] * g_max_by_z_tech[z,ndisp] for z in Z, ndisp in NONDISP)

feed_in = Dict(z => sum(feed_in_by_z_nondisp[z,ndisp] for ndisp in NONDISP, dim=2) for z in Z)

### other utility ###
successor(arr, x) = x == length(arr) ? 1 : x + 1

### ntc ###
ntc = dictzip(lines, [:from, :to] => :line_capacity)
for (k,v) in ntc ntc[reverse(k)] = v end

# actual model creation
###############################################################################
m = Model(Clp.Optimizer)

@variables m begin
    g_max[disp] >= G[disp=P_DISP, T] >= 0
    CU[Z,T] >= 0
    g_max[s] >= D_stor[s=P_S,T] >= 0
    stor_max[s] >= L_stor[s=P_S,T] >= 0
    H[P_CHP, T] >= 0
    get(ntc, (z,zz), 0) >= FLOW[z=Z,zz=Z,T] >= 0
end

@objective(m, Min,
    sum(mc[map_id2tech[disp]] * G[disp,t] for disp in P_DISP, t in T) );
@constraint(m, ElectricityBalance[z=Z, t=T],
    sum(G[disp,t] for disp in intersect(map_zone2plist[z],P_DISP))
    + feed_in[z][t]
    - sum(D_stor[s,t] for s in intersect(map_zone2plist[z],P_S))
    + sum(FLOW[zz,z,t]-FLOW[z,zz,t] for zz in Z)
    - CU[z,t]
    ==
    elec_demand[z][t])

@constraint(m, HeatBalance[z=Z,t=T],
    sum(H[ht, t] for ht in intersect(map_zone2plist[z],P_CHP))
    >=
    heat_demand[z][t])

@constraint(m, CoGeneration[chp=P_CHP, t=T],
    G[chp,t] == heat_ratio[map_id2tech[chp]] * H[chp, t])

@constraint(m, StorageLevel[s=P_S, t=T],
    L_stor[s, successor(T,t)]
    ==
    L_stor[s, t]
    + sqrt(eff[map_id2tech[s]])*D_stor[s,t]
    - (1/sqrt(eff[map_id2tech[s]])) * G[s,t])

optimize!(m)

# post processing
#############################################################################
generic_names(len::Int) = [Symbol("x$i") for i in 1:len]
function DataFrame(x::JuMP.Containers.DenseAxisArray,
    names::AbstractVector{Symbol}=generic_names(length(x.axes));
    kwargs...)

    if length(names) != length(x.axes)
        throw(ArgumentError("Length of argument 'names' is $(length(names)) and
        does not fit the variable's dimension of $(length(x.axes))"))
    end
    push!(names, :value)

    iter = Iterators.product(x.axes...) |> collect |> vec
    df = DataFrame([(i..., value(x[i...])) for i in iter])
    rename!(df, names)
    return df
end

result_G = DataFrame(G, [:id, :hour])
insertcols!(result_G, 2, :zone => [map_id2z[id] for id in result_G[!,:id]])
insertcols!(result_G, 3, :technology => [map_id2tech[id] for id in result_G[!,:id]])

gen_by_tech = combine(groupby(result_G, [:zone, :technology, :hour]),
    :value => sum => :value)

result_feedin = DataFrame(
        (id=id,
        zone=map_id2z[id],
        technology=map_id2tech[id],
        hour=t,
        value=feed_in_by_z_nondisp[map_id2z[id],map_id2tech[id]][t]
        )
    for id in P_NONDISP, t in T)

gen_by_tech2 = combine(groupby(result_feedin, [:zone, :technology, :hour]),
    :value => sum => :value)

result_CU = DataFrame(CU, [:zone, :hour])
result_CU[!,:technology] .= "curtailment"

result_D = DataFrame(D_stor, [:id, :hour])
insertcols!(result_D, 2, :zone => [map_id2z[id] for id in result_D[!,:id]])
insertcols!(result_D, 3, :technology => [map_id2tech[id] for id in result_D[!,:id]])

demand_by_storage = combine(groupby(result_D, [:zone, :technology, :hour]),
    :value => sum => :value)
df_elec_demand = DataFrame((zone=k, technology="demand", hour=x, value=v[x])
    for (k,v) in elec_demand for x in 1:length(v))

result_generation = vcat(gen_by_tech2, gen_by_tech)
result_demand = vcat(demand_by_storage, result_CU, df_elec_demand)


### Plotting
colordict = Dict("pv" => :yellow, "wind" => :lightblue, "pumped_hydro" => :darkblue,
    "battery" => :lightgrey, "p1" => :brown, "p2" => :grey, "demand" => :darkgrey,
    "curtailment" => :red)

function plot_energybalance(df_gen::DataFrame, df_dem::DataFrame, z::AbstractString)

    df_generation = filter(x-> x.zone == z, df_gen)
    df_demand = filter(x-> x.zone == z, df_dem)

    table_gen = unstack(df_generation, :hour, :technology, :value)
    labels = names(table_gen[:,Not(:hour)]) |> permutedims
    colors = [colordict[tech] for tech in labels]
    data_gen = Array(table_gen[:,Not(:hour)])

    p = areaplot(data_gen, label=labels, color=colors, width=0,
        leg=:outertopright, title="Dispatch of zone $z")

    table_dem = unstack(df_demand, :hour, :technology, :value)
    labels = names(table_dem[:,Not(:hour)]) |> permutedims
    colors = [colordict[tech] for tech in labels]
    data_dem = - Array(table_dem[:,Not(:hour)])

    areaplot!(p, data_dem, label=labels, color=colors, width=0)
    hline!(p, [0], color=:black, label="", width=2)

    return p
end

z1 = plot_energybalance(result_generation, result_demand, "z1")
#savefig("results_z1.pdf")
z2 = plot_energybalance(result_generation, result_demand, "z2")
#savefig("results_z2.pdf")
