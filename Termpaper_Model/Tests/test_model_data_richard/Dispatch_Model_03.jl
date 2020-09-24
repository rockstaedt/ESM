using JuMP
using Gurobi
using DataFrames
using LinearAlgebra
using CSV
using Plots
using Statistics

###############################################################################
### Data Import
###############################################################################

data_path = "data"


P_df = joinpath(data_path,"plants_test.csv") |> CSV.read
D_df = joinpath(data_path,"demand_test.csv") |> CSV.read
R_df = joinpath(data_path,"avail_test.csv") |> CSV.read;


#to sum up demand
D_df = select!(D_df, Not(:Column1))
D_df = select!(D_df, Not(:index))

###############################################################################
### Data Preprocessing
###############################################################################

Demand = sum(eachcol(D_df));
# dispatchable plants
DISP = P_df[P_df.plant_type .== "conventional",:].index
append!(DISP, P_df[P_df.plant_type .== "hydro_ror",:].index)

# nondispatchable plants
# excluding "other_res" because no availability data
NDISP = P_df[P_df.plant_type .== "wind onshore",:].index
append!(NDISP, P_df[P_df.plant_type .== "wind offshore",:].index)
append!(NDISP, P_df[P_df.plant_type .== "solar",:].index)

PLANTS = vec(Array(select(P_df,"index")))

T = T = 1:size(D_df, 1) |> collect

DISP_PLANTS = copy(P_df)
NONDISP_PLANTS = copy(P_df)
delete!(DISP_PLANTS,1:size(P_df,1))
delete!(NONDISP_PLANTS,1:size(P_df,1))

for i in DISP
    f = P_df[P_df.index .==i,:]
    append!(DISP_PLANTS,f)
end
for i in NDISP
    f = P_df[P_df.index .==i,:]
    append!(NONDISP_PLANTS,f)
end
DISP_PLANTS
NONDISP_PLANTS;

# map marginal costs to conventional power plants
mc_c = Dict(row.index => row.mc_el for row in eachrow(DISP_PLANTS))
# map max of generation of conventional power plants
gmax_c = Dict(row.index => row.g_max for row in eachrow(DISP_PLANTS));

res_feed_in = Dict()
for row in eachrow(NONDISP_PLANTS)
        for t in T
            res_feed_in[row.index,t] = row.g_max * R_df[t,row.index]
        end
end
res_feed_in





####MODEL

Dax = Model(Gurobi.Optimizer)
@variables Dax begin
    G[DISP,T] >= 0
    Gr[NDISP,T] >= 0
    CU[T]>=0
end

@objective(Dax, Min, sum(mc_c[disp] * G[disp,timestep] for disp in DISP, timestep in T));


@constraint(Dax, Max_Generation[plant=DISP,timestep = T],
            G[plant,timestep] <= gmax_c[plant]);

@constraint(Dax, Max_Generation_Res[plant=NDISP,timestep = T],
            Gr[plant,timestep] == res_feed_in[plant,timestep]);

@constraint(Dax, EnergyBalance[timestep = T],
            sum(G[plant,timestep] for plant in DISP)
            + sum(Gr[ndisp,timestep] for ndisp in NDISP)
            - CU[timestep] == Demand[timestep]);


JuMP.optimize!(Dax)
JuMP.objective_value(Dax)


Price=dual.(EnergyBalance).data
mean(Price)
maximum(Price)
minimum(Price)
maximum(value.(CU))
