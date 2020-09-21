using JuMP
using Gurobi
using DataFrames
using LinearAlgebra
using CSV
using Plots

###############################################################################
### Data Import
###############################################################################

data_path = "data"
P_df = joinpath(data_path,"plants_test.csv") |> CSV.read
D_df = joinpath(data_path,"demand_test.csv") |> CSV.read
R_df = joinpath(data_path,"avail_test.csv") |> CSV.read;

#Demand
D_df = select!(D_df, Not(:Column1))
Demand = sum(eachcol(D_df))     #demand array
Demand_df = DataFrame(A = Demand) #demand dataframe


###Split of res plants and conventional plants###
Res_Array = []
Disp_Array = []
mc = Dict()
gmax_c = Dict()
gmax_r = Dict()
for i in eachrow(P_df)

    if i.fuel =="sun"
        push!( Res_Array, i.index  )
        gmax_r[i.index]=i.g_max
    elseif i.fuel =="wind"
        push!( Res_Array, i.index  )
        gmax_r[i.index]=i.g_max
    else
        push!( Disp_Array, i.index  )
        mc[i.index]=i.mc_el
        gmax_c[i.index]=i.g_max
    end
end
Disp_Array;
mc;
gmax_r;
T = collect(1:nrow(D_df));

###Res Dict
Rdict = Dict()
for r in Res_Array
    for t in T
        Rdict[r,t] = R_df[t,r]*gmax_r[r]


    end
end
Rdict;


###Model
Dax = Model(Gurobi.Optimizer)
@variables Dax begin
    G[Disp_Array,T] >= 0
    Gr[Res_Array,T] >= 0
end

@objective(Dax, Min, sum(mc[disp] * G[disp,timestep] for disp in Disp_Array, timestep in T));


@constraint(Dax, Max_Generation[plant=Disp_Array,timestep = T],
            G[plant,timestep] <= gmax_c[plant]);

@constraint(Dax, Max_Generation_Res[plant=Res_Array,timestep = T],
            Gr[plant,timestep] <= Rdict[plant,timestep]);

@constraint(Dax, EnergyBalance[timestep = T],
            sum(G[plant,timestep] for plant in Disp_Array)
            + sum(Gr[ndisp,timestep] for ndisp=Res_Array) == Demand_df[T,1]);


JuMP.optimize!(Dax)
JuMP.objective_value(Dax)
