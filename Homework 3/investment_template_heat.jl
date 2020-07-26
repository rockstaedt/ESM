using JuMP
using Clp
using Plots
using DataFrames, CSV

data_path = "data"
time_series = CSV.read(joinpath(data_path, "timedata.csv"))
# add heat time series
insertcols!(time_series, :heat => CSV.read(joinpath(data_path, "timedata_with_heat.csv"))[!,:heat]|> Vector)
tech_data = CSV.read(joinpath(data_path, "technologies.csv"))
# read new tech
tech_new_data = CSV.read(joinpath(data_path, "new_tech.csv"))
# add column dispatchable
insertcols!(tech_new_data, 2, :dispatchable => 1)
# add new tech to tech data frame as row
tech_data = vcat(tech_data, tech_new_data)
tech_data



### data preprocessing ###
T = 1:size(time_series, 1) |> collect
P = tech_data[:,:technology] |> Vector
DISP = tech_data[tech_data[!,:dispatchable] .== 1 ,:technology]
NONDISP = tech_data[tech_data[!,:dispatchable] .== 0 ,:technology]
S = tech_data[tech_data[!,:investment_storage] .> 0 ,:technology]

### parameters ###
annuity_factor(n,r) = r * (1+r)^n / (((1+r)^n)-1)

interest_rate = 0.04
ic_generation_cap = Dict{String, Float64}()
ic_charging_cap = Dict{String, Float64}()
ic_storage_cap = Dict{String, Float64}()
eff_in = Dict{String, Float64}()
eff_out = Dict{String, Float64}()
vc = Dict{String, Float64}()


for row in eachrow(tech_data)
    af = annuity_factor(row.lifetime, interest_rate)
    ic_generation_cap[row.technology] = row.investment_generation * af

    iccc = row.investment_charge * af
    iccc > 0 && (ic_charging_cap[row.technology] = iccc)

    icsc = row.investment_storage * af
    icsc > 0 && (ic_storage_cap[row.technology] = icsc)

    row.storage_efficiency_in > 0 && (eff_in[row.technology] = row.storage_efficiency_in)
    row.storage_efficiency_out > 0 && (eff_out[row.technology] = row.storage_efficiency_out)

    vc[row.technology] = row.vc
end

demand = time_series[:,:demand] |> Array

### HEAT SECTION

# add heat demand
heat_demand = time_series[:,:heat] |> Array

availability = Dict(nondisp => time_series[:,nondisp] for nondisp in NONDISP)

successor(arr, x) = x == length(arr) ? 1 : x + 1

dispatch_scale = 8760/length(T)

### model ###
m = Model(Clp.Optimizer)

@variables m begin
    # variables from our dispatch model
    G[DISP, T] >= 0
    CU[T] >= 0
    D_stor[S,T] >= 0
    L_stor[S,T] >= 0

    H[P2H, T] >= 0
    H_P2H[P2H, T] >= 0

    # new variables for our investment model
    CAP_G[P] >= 0
    CAP_D[S] >= 0
    CAP_L[S] >= 0
end

@objective(m, Min,
    sum(vc[disp] * G[disp,t] for disp in DISP, t in T) * dispatch_scale
    + sum(ic_generation_cap[p] * CAP_G[p] for p in P)
    + sum(ic_charging_cap[s] * CAP_D[s] for s in S if haskey(ic_charging_cap, s))
    + sum(ic_storage_cap[s] * CAP_L[s] for s in S)
)

@expression(m, feed_in[ndisp=NONDISP, t=T], availability[ndisp][t]*CAP_G[ndisp])

@constraint(m, ElectricityBalance[t=T],
    sum(G[disp,t] for disp in DISP)
    + sum(feed_in[ndisp,t] for ndisp in NONDISP)
    - sum(D_stor[s,t] for s in S)
    - sum(H_P2H[p2h,t] for p2h in P2H)
    - CU[t]
    ==
    demand[t])

@constraint(m, MaxGeneration[disp=DISP, t=T],
    G[disp,t] <= CAP_G[disp])

@constraint(m, MaxCharge[s=S, t=T; haskey(ic_charging_cap, s)],
    D_stor[s,t] <= CAP_D[s])

@constraint(m, SymmetricChargingPower[s=S, t=T; !(haskey(ic_charging_cap, s))],
    CAP_G[s] == CAP_D[s])

@constraint(m, MaxLevel[s=S, t=T],
    L_stor[s,t] <= CAP_L[s])

@constraint(m, StorageLevel[s=S, t=T],
    L_stor[s, successor(T,t)]
    ==
    L_stor[s, t]
    + eff_in[s]*D_stor[s,t]
    - (1/eff_out[s]) * G[s,t] )

@constraint(m, HeatBalance[t=T],
    sum(H[ht, t] for ht in HT)
    >=
    heat_demand[t])

@constraint(m, Power2Heat[p2h=P2H, t=T],
    H[p2h,t] == cop[p2h] * H_P2H[p2h, t])

@constraint(m, MaxHeat[disp=DISP, t=T],
    H[disp,t] <= CAP_G[disp])

optimize!(m)


######
generic_names(len::Int) = [Symbol("x$i") for i in 1:len]
function get_result(x::JuMP.Containers.DenseAxisArray,
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

colordict = Dict("pv" => :yellow, "wind" => :lightblue, "seasonal_storage" => :darkblue,
    "battery" => :lightgrey, "demand" => :darkgrey, "curtailment" => :brown,
    "power-to-heat" => :red)

######## plot electricity balance ###########

result_G = get_result(G, [:technology, :hour])
result_feed_in = get_result(feed_in, [:technology, :hour])

result_charging = get_result(D_stor, [:technology, :hour])
result_CU = get_result(CU, [:hour])
result_CU[!,:technology] .= "curtailment"
df_demand = DataFrame(hour=T, technology="demand", value=demand)


result_generation = vcat(result_feed_in, result_G)
result_demand = vcat(result_charging, result_CU, df_demand)

table_gen = unstack(result_generation, :hour, :technology, :value)
table_gen = table_gen[!,[NONDISP..., DISP...]]
labels = names(table_gen) |> permutedims
colors = [colordict[tech] for tech in labels]
data_gen = Array(table_gen)

balance_plot = areaplot(data_gen, label=labels, color=colors, width=0,
    leg=:outertopright)

table_dem = unstack(result_demand, :hour, :technology, :value)
# add power-to-heat accordingly here
#table_dem = table_dem[!,["demand","power-to-heat", S...,"curtailment"]]
table_dem = table_dem[!,["demand", S...,"curtailment"]]
labels2 = names(table_dem) |> permutedims
colors2 = [colordict[tech] for tech in labels2]
replace!(labels2, [item => "" for item in intersect(labels2, labels)]...)
data_dem = -Array(table_dem)

areaplot!(balance_plot, data_dem, label=labels2, color=colors2, width=0,
    leg=:outertopright)

hline!(balance_plot, [0], color=:black, label="", width=2)

#################################

df_installed_gen = get_result(CAP_G, [:technology])
x = df_installed_gen[!,:technology]
y = df_installed_gen[!,:value]
c = [colordict[tech] for tech in x]
p1 = bar(x, y, color=c, leg=false, title="Installed power generation",
    ylabel="MW", rotation=90)

df_installed_charge = get_result(CAP_D, [:technology])
x = df_installed_charge[!,:technology]
y = df_installed_charge[!,:value]
c = [colordict[tech] for tech in x]
p2 = bar(x, y, color=c, leg=false, title="Installed power charging", rotation=90,
    ylim=ylims(p1))

df_installed_storage = get_result(CAP_L, [:technology])
x = df_installed_storage[!,:technology]
y = df_installed_storage[!,:value]
c = [colordict[tech] for tech in x]
p3 = bar(x, y, color=c, leg=false, title="Installed storage capacity",
    ylabel="MWh", rotation=90)

plot(p1,p2,p3,layout=(1,3), titlefontsize=6, tickfontsize=6)
