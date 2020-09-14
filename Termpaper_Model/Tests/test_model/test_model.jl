###############################################################################
### Packages
###############################################################################

using JuMP
using Gurobi
using DataFrames
using LinearAlgebra
using CSV
using Plots

## hallo Theo

# Script to calculate PTDF Matrix
include("ptdf_calculation.jl")

###############################################################################
### Data Import
###############################################################################

data_path = "data"

L_df = joinpath(data_path,"lines.csv") |> CSV.read
N_df = joinpath(data_path,"nodes.csv") |> CSV.read
P_df = joinpath(data_path,"plants.csv") |> CSV.read
D_df = joinpath(data_path,"demand.csv") |> CSV.read
R_df = joinpath(data_path,"res.csv") |> CSV.read
S_df = joinpath(data_path,"storages.csv") |> CSV.read

###############################################################################
### Data Preprocessing
###############################################################################

## PTDF Matrix

b = vec(Array(select(L_df,"b")))
A = zeros(nrow(L_df),nrow(N_df))
NODES = vec(Array(select(N_df,"index")))
SLACK = 14

# calculation of incidence matrix A
counter= 1
for i in eachrow(L_df)
    wert1 = L_df[counter,:][2]
    wert2 = L_df[counter,:][3]
    A[counter,wert1] =1
    A[counter,wert2] =-1
    global counter +=1
end
A

ptdf = calculate_ptdf(A,b,NODES,SLACK)
ptdf

# dispatchable plants
DISP = P_df[P_df.renewable .=="no",:].index
# nondispatchable plants, shouldnt be that important
NDISP = P_df[P_df.renewable .=="yes",:].index
PLANTS = vec(Array(select(P_df,"index")))
LINES = vec(Array(select(L_df,"index"))) #lines
l_max = vec(Array(select(L_df,"lmax"))) #line maxes
T = [1,2] #timesteps
demand= D_df #demand dataframe


#df for disp plants

#make empty copy
DISP_PLANTS = copy(P_df)
NONDISP_PLANTS = copy(P_df)
deleterows!(DISP_PLANTS,1:16)
deleterows!(NONDISP_PLANTS,1:16)

for i in DISP
    f = P_df[P_df.index .==i,:]
    append!(DISP_PLANTS,f)
end
for i in NDISP
    f = P_df[P_df.index .==i,:]
    append!(NONDISP_PLANTS,f)
end
DISP_PLANTS
NONDISP_PLANTS

plants_at_node = Dict()
mc = Dict()
gmax_c = Dict()
for i in eachrow(P_df)
    if !haskey(plants_at_node,i["node"])
        plants_at_node[i["node"]]=[i["index"]]
    else
        push!(plants_at_node[i["node"]],i["index"])
    end
    mc[i["index"]] = i["mc"]
    gmax_c[i["index"]] = i["g_max_c"]
end
plants_at_node[SLACK]=[]
plants_at_node

#to enter slack
Kopie = copy(P_df)
push!(Kopie,["0" 14 "0" 0 0 "0" 0])

#gmax_r für alle nodes, timesteps
#Capable for two res techs in a single node(pv and wind must be aggregated before hand)
gmax_r = Dict()
Techs = [1,2,3]#1=wind,2=pv,3=other



for n in NODES
    for t in T
        for i in Techs
            gmax_r[n,t,i] = 0
        end
        for i in 1:nrow(P_df[P_df.node .==n,:])
            if P_df[P_df.node .==n,:][i,:]["tech"] == "wind"
                gmax_r[n,t,1]= P_df[P_df.node .==n,:][i,:]["g_max_r"]*R_df[t,2]
            elseif P_df[P_df.node .==n,:][i,:]["tech"] == "pv"
                gmax_r[n,t,2] = P_df[P_df.node .==n,:][i,:]["g_max_r"]*R_df[t,3]
            end
        end
    end
end

NODES = vec(Array(select(N_df,"index")))
NODESCOPY = copy(NODES)
for i in S_df.node
    filter!(x-> x !=i ,NODESCOPY)
end

#Storage Data Preparation
storage_at_node = Dict()
gmax_s = Dict()
mcs=Dict()
for i in eachrow(S_df)
    if !haskey(storage_at_node,i["node"])
        storage_at_node[i["node"]]=[i["index"]]
    else
        push!(storage_at_node[i["node"]],i["index"])
    end
    gmax_s[i["index"]] = i["g_max_s"]
    mcs[i["index"]] = i["mc"]

end

storage_at_node
STORAGES = vec(Array(select(S_df,"index")));
gmax_s;#storage cap and g max are equal in each timestep
mcs;

for n in NODESCOPY
    storage_at_node[n] =[]
end

T=[1,2]
Ts = copy(T)
push!(Ts,length(T)+1)# Timesteps for storage levels


###############################################################################
### Model
###############################################################################

Fax = Model(Gurobi.Optimizer)
@variables Fax begin
    G[PLANTS,T] >= 0
    INJ[NODES,T]
    Gr[NODES,T] >= 0
    D[STORAGES,Ts]>=0
    Gs[STORAGES,Ts]>=0
    L[STORAGES,Ts]>=0

end

@objective(Fax, Min, sum(mc[disp] * G[disp,timestep] for disp in DISP, timestep in T)
                    +sum(mcs[s]*Gs[s,timestep] for s in STORAGES,timestep in T ));

@constraint(Fax, Max_Generation[plant=PLANTS,timestep = T],
            G[plant,timestep] <= gmax_c[plant]);

@constraint(Fax, Max_Gen_Res[node=NODES,timestep=T],
            Gr[node,timestep] <= sum(gmax_r[node,timestep,tech] for tech in Techs));
#so far so good



@constraint(Fax, EnergyBalance[node=NODES,timestep = T],
            sum(G[plant,timestep] for plant in plants_at_node[node])
            - demand[timestep,node]  + Gr[node,timestep]
            - sum(D[s,timestep] for s in storage_at_node[node])
            + sum(Gs[s,timestep] for s in storage_at_node[node])== INJ[node,timestep]);


#storage constraints

@constraint(Fax,Input_Cap[s=STORAGES,timestep=T],
            D[s,timestep] <= gmax_s[s]);

@constraint(Fax,Generation_Cap[s=STORAGES,timestep=T],
            Gs[s,timestep] <= gmax_s[s]);

@constraint(Fax,Storage_Cap[s=STORAGES,timestep=T],
            L[s,timestep] <= gmax_s[s]);

@constraint(Fax,Storage_Level[s=STORAGES,timestep=2:length(T)],
            L[s,timestep+1]==L[s,timestep]-Gs[s,timestep]+D[s,timestep]);

@constraint(Fax,Start_Level[s=STORAGES],
            L[s,1]==gmax_s[s]);


#Line Constraints
@constraint(Fax, LineMax[line=LINES],
            sum(ptdf[line, node] * INJ[node,timestep] for node in NODES,timestep in T) <= l_max[line]
            );

@constraint(Fax, LineMin[line=LINES],
            sum(ptdf[line, node] * INJ[node,timestep] for node in NODES,timestep in T) >= -l_max[line]
            );

@constraint(Fax, sum(INJ[node,timestep] for node in NODES,timestep in T) == 0);

JuMP.optimize!(Fax)
JuMP.objective_value(Fax)


###############################################################################
### Results
###############################################################################

export_path = "export_files"

generation_t1 = DataFrame(Plants = PLANTS,
                        t = "t1",
                        DISP = value.(G).data[:,1]
                        )
generation_t2 = DataFrame(Plants = PLANTS,
                        t = "t2",
                        DISP = value.(G).data[:,2]
                        )

generation = vcat(generation_t1, generation_t2)


CSV.write(joinpath(export_path, "generation.csv"), generation)

value.(Gr)

value.(Gs)

value.(D)

value.(L)

Price=dual.(EnergyBalance).data
