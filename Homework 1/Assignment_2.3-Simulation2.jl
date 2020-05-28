using JuMP
using Clp
using Plots
using DataFrames, CSV

gr()

Plots.GRBackend()

data = CSV.read("timedata.csv")
rename!(data, :solar => :pv)

T = Vector(data[:,:hour])
P = ["p1","p2","pv","wind"]
DISP = ["p1","p2"]
NONDISP = ["pv","wind"]

demand = Dict(data[:,:hour] .=> data[:,:demand]./1000)

mc = Dict("p1" =>  10, "p2" => 20)
g_max = Dict("p1" =>  25, "p2" => 40)

pv_installed = 50
wind_installed = 60

res_feed_in = Dict((res, row["hour"]) => row[res] * cap
    for (res,cap) in [("pv",pv_installed),("wind",wind_installed)], row in eachrow(data))

S = ["PumpedHydro", "Battery"]

eff = Dict("PumpedHydro" => 1, "Battery" => 1)

mc["PumpedHydro"] = 1
mc["Battery"] = 1

g_max["PumpedHydro"] = 8
g_max["Battery"] = 2

stor_max = Dict("PumpedHydro" => 50, "Battery" => 4)

P = vcat(P, S)

DISP = vcat(DISP, S)

m = Model(Clp.Optimizer)

@variables m begin
    G[DISP, T] >= 0
    CU[T] >= 0
    D_stor[S,T] >= 0
    L_stor[S,T] >= 0
end

@objective(m, Min, sum(mc[disp] * G[disp,t] for disp in DISP, t in T))

@constraint(m, EnergyBalance[t=T],
    sum(G[disp,t] for disp in DISP)
    + sum(res_feed_in[ndisp,t] for ndisp=NONDISP)
    - sum(D_stor[s,t] for s in S)
    - CU[t]
    ==
    demand[t])

@constraint(m, MaxGeneration[disp=DISP, t=T],
    G[disp,t] <= g_max[disp])

@constraint(m, MaxCharge[s=S, t=T],
    D_stor[s,t] <= g_max[s])

@constraint(m, StorageLevel[s=S, t=1:length(T)-1],
    L_stor[s, T[t + 1]]
    ==
    L_stor[s, T[t]]
    + sqrt(eff[s])*D_stor[s,T[t]]
    - (1/sqrt(eff[s])) * G[s,T[t]])

last_t = T[end]
@constraint(m, StorageLevelEnd[s=S],
    L_stor[s, T[1]]
    ==
    L_stor[s, last_t]
    + sqrt(eff[s])*D_stor[s,last_t]
    - (1/sqrt(eff[s])) * G[s,last_t])

optimize!(m)

total_cost = round(objective_value(m), digits=2)

result_G = value.(G)

feedin = [res_feed_in[ndisp,t] for ndisp in NONDISP, t in T]

generation = vcat(feedin, result_G)

curtailment = value.(CU)

d = [demand[t] for t in T]

areaplot(
    generation',
    label=["PV" "Wind" "P1" "P2" "PumpedHydro" "Battery"],
    color=[:yellow :lightblue :brown :grey :blue :purple],
    legend=:bottomright,
    xlabel="Hours",
    ylabel="GW",
    width=0,
    title="Simulation 2 - Total Cost= $total_cost")

plot!(d, color=:black, width=3, label = "Demand")
plot!(curtailment.data, color=:black, width=2, label="Curtailment", linestyle=:dash)

stor = -value.(D_stor)

areaplot!(
    stor.data',
    label="",
    color=[:blue :purple],
    width=0)

hline!([0], color=:black, width=2, label="")

price = dual.(EnergyBalance)

plot!(twinx(), price.data, color=:red, width=2, leg=false, ylabel="Cost per GW")

# savefig("figure-simulation2.pdf")
