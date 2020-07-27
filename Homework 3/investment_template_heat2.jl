using JuMP
using Clp
using Plots
using DataFrames, CSV

data_path = "data"
time_series = CSV.read(joinpath(data_path, "timedata.csv"))
# add heat time series
insertcols!(time_series, :heat => CSV.read(joinpath(data_path, "timedata_with_heat.csv"))[!,:heat]|> Vector)

tech_data = CSV.read(joinpath(data_path, "technologies.csv"))
# read heat tech data
heat_tech_data = CSV.read(joinpath(data_path, "new_tech.csv"))



### data preprocessing ###
T = 1:size(time_series, 1) |> collect
P = tech_data[:,:technology] |> Vector
DISP = tech_data[tech_data[!,:dispatchable] .== 1 ,:technology]
NONDISP = tech_data[tech_data[!,:dispatchable] .== 0 ,:technology]
S = tech_data[tech_data[!,:investment_storage] .> 0 ,:technology]
H = ["p2h"]

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

# add investment variables for heating techs

ic_generation_cap_h = Dict{String, Float64}()
ic_storage_cap_h = Dict{String, Float64}()
eff_in_h = Dict{String, Float64}()
eff_out_h = Dict{String, Float64}()
vc_h = Dict{String, Float64}()

af_h = annuity_factor(heat_tech_data.lifetime[1], interest_rate)
ic_generation_cap_h[heat_tech_data.technology[1]] = heat_tech_data.investment_generation[1] * af_h

icsc = heat_tech_data.investment_storage[1] * af_h
icsc > 0 && (ic_storage_cap_h[heat_tech_data.technology[1]] = icsc)

heat_tech_data.storage_efficiency_in[1] > 0 && (eff_in_h[heat_tech_data.technology[1]] = heat_tech_data.storage_efficiency_in[1])
heat_tech_data.storage_efficiency_out[1] > 0 && (eff_out_h[heat_tech_data.technology[1]] = heat_tech_data.storage_efficiency_out[1])

vc_h[heat_tech_data.technology[1]] = heat_tech_data.vc[1]

# demand
demand = time_series[:,:demand] |> Array

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
    # heat
    H_G[H, T] >= 0
    H_P2H[H,T] >= 0
    H_D_S[H,T] >= 0
    H_L_S[H, T] >= 0

    # new variables for our investment model
    CAP_G[P] >= 0
    CAP_D[S] >= 0
    CAP_L[S] >= 0
    # heat
    CAP_H_G[H] >= 0
    CAP_H_D[H] >= 0
    CAP_H_L[H] >= 0
end

@objective(m, Min,
    sum(vc[disp] * G[disp,t] for disp in DISP, t in T) * dispatch_scale
    + sum(ic_generation_cap[p] * CAP_G[p] for p in P)
    + sum(ic_charging_cap[s] * CAP_D[s] for s in S if haskey(ic_charging_cap, s))
    + sum(ic_storage_cap[s] * CAP_L[s] for s in S)
    + sum(vc_h[h] * H_G[h,t] for h in H, t in T) * dispatch_scale
    + sum(ic_generation_cap_h[h] * CAP_H_G[h] for h in H)
    + sum(ic_storage_cap_h[h] * CAP_H_L[h] for h in H)
)

@expression(m, feed_in[ndisp=NONDISP, t=T], availability[ndisp][t]*CAP_G[ndisp])

@constraint(m, ElectricityBalance[t=T],
    sum(G[disp,t] for disp in DISP)
    + sum(feed_in[ndisp,t] for ndisp in NONDISP)
    - sum(D_stor[s,t] for s in S)
    - sum(H_P2H[h, t] for h in H)
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

@constraint(m, MaxGenerationHeat[h=H, t=T],
    H_G[h,t] <= CAP_H_G[h])

@constraint(m, MaxLevelHeat[h=H, t=T],
    H_L_S[h,t] <= CAP_H_L[h])

@constraint(m, StorageLevel[s=S, t=T],
    L_stor[s, successor(T,t)]
    ==
    L_stor[s, t]
    + eff_in[s]*D_stor[s,t]
    - (1/eff_out[s]) * G[s,t])

@constraint(m, StorageLevelHeat[h=H, t=T],
    H_L_S[h, successor(T,t)]
    ==
    H_L_S[h, t]
    + eff_in_h[h]*H_D_S[h,t]
    - (1/eff_out_h[h]) * H_G[h,t])

@constraint(m, HeatBalance[t=T],
    sum(H_G[h, t] for h in H)
    >=
    heat_demand[t])

@constraint(m, Power2Heat[h=H, t=T],
    H_G[h,t] == eff_in_h[h] * H_P2H[h, t])

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
    "p2h" => :red)

######## plot electricity balance ###########

result_G = get_result(G, [:technology, :hour])
result_feed_in = get_result(feed_in, [:technology, :hour])

result_charging = get_result(D_stor, [:technology, :hour])
result_CU = get_result(CU, [:hour])
result_CU[!,:technology] .= "curtailment"
df_demand = DataFrame(hour=T, technology="demand", value=demand)
df_heat_g = get_result(H_P2H, [:technology, :hour])
#df_heat_d_s = get_result(H_D_S, [:technology, :hour])


result_generation = vcat(result_feed_in, result_G)
result_demand = vcat(result_charging, result_CU, df_demand, df_heat_g)

table_gen = unstack(result_generation, :hour, :technology, :value)
table_gen = table_gen[!,[NONDISP..., DISP...]]
labels = names(table_gen) |> permutedims
colors = [colordict[tech] for tech in labels]
data_gen = Array(table_gen)

balance_plot = areaplot(data_gen, label=labels, color=colors, width=0,
    leg=:outertopright)

table_dem = unstack(result_demand, :hour, :technology, :value)
# add power-to-heat accordingly here
table_dem = table_dem[!,["demand","p2h", S...,"curtailment"]]
#table_dem = table_dem[!,["demand", S...,"curtailment"]]
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

####
heat_generation = value.(H_L_S).data'

heat_balance = areaplot(heat_generation,
    color=[:darkgrey :darkred :orange],
    width=0,
    legend=false,
    yaxis="GW")

plot!(heat_demand, color=:black, width=2, label="")
