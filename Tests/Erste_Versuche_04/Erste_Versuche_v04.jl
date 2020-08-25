using JuMP
using Clp
using DataFrames
using LinearAlgebra
using CSV
L_df = CSV.read("lines.csv")
N_df = CSV.read("nodes.csv")
P_df = CSV.read("plants.csv")
D_df = CSV.read("demand.csv")

b = vec(Array(select(L_df,"b")))
A = zeros(nrow(L_df),nrow(N_df))
NODES = vec(Array(select(N_df,"index")))
SLACK = 14
PLANTS = vec(Array(select(P_df,"index")))
LINES = vec(Array(select(L_df,"index")))
T = [1,2]

demand= D_df
#Demand = Dict()
#Demand["t1"]=Array(D_df[1,:])
#Demand["t2"]=Array(D_df[2,:])
l_max = vec(Array(select(L_df,"lmax")))


#calculation of incidence matrix A
counter= 1
for i in eachrow(L_df)
    wert1 = L_df[counter,:][2]
    wert2 = L_df[counter,:][3]
    A[counter,wert1] =1
    A[counter,wert2] =-1
    global counter +=1
end
A

#Account for double entries
#plants_at_node dict
#gmax dict
#mc
plants_at_node = Dict()
gmax = Dict()
mc = Dict()

for i in eachrow(P_df)
    if !haskey(plants_at_node,i["node"])
        plants_at_node[i["node"]]=[i["index"]]
    else
        push!(plants_at_node[i["node"]],i["index"])
    end
    gmax[i["index"]] = i["g_max"]
    mc[i["index"]] = i["mc"]



end
plants_at_node[SLACK]=[]


#calculate ptdf matrix source:lecture
function calculate_ptdf(A, b_vector, NODES, SLACK)
    """ Create PTDF Based on:
    Topology (A, incedence Matrix): LxN Matrix (1, -1) for nodes connected by line
    Susceptance Vector (B): Linesusceptance, calculated by line paramters and propotional to length
    NODES: Array of nodes
    SLACK: Reference Node to create the inversse of A.T * Bd A

    From documentation: PTDF = (Bd*A)(A.T * Bd * A)^-1
    Bd: Suceptance Vector (for each line) on Diagonal Matrix (LxL)
    Bl: Line Susceptance Matrix Bd*A
    Bn: Nodes Susceptance Matrix (A'*B)*A
    B_inv = Bn^-1 (without slack, otherwise singular)
    PTDF = Bl * B_inv
    """
    ## Calculation based on indices rather than Node Names (of type Str or Symb)
    idx_nodes = collect(1:length(NODES))
    idx_slack = findfirst(x -> x==SLACK, NODES)
    #A = incedence
    Bd = Diagonal(b_vector)
    Bl = Bd*A # Line suceptance Matrix
    Bn = (A'*Bd)*A # Nodes suceptance Matrix

    B_inv = zeros(length(NODES),length(NODES))
    B_inv[setdiff(idx_nodes,idx_slack), setdiff(idx_nodes,idx_slack)] = inv(Bn[setdiff(idx_nodes,idx_slack), setdiff(idx_nodes,idx_slack)])
    PTDF = Bl*B_inv
    return PTDF
end


ptdf = calculate_ptdf(A,b,NODES,SLACK)

#ptdf_dict = Dict((l, n) => ptdf[l_idx, n_idx] for (n_idx, n) in enumerate(NODES),
 #                                                 (l_idx, l) in enumerate(LINES))




 ##Model and stuff
Fax = Model(Clp.Optimizer)
@variables Fax begin
    G[PLANTS,T] >= 0
    INJ[NODES,T]
end

@objective(Fax, Min, sum(mc[plant] * G[plant,timestep] for plant in PLANTS, timestep in T));

@constraint(Fax, Max_Generation[plant=PLANTS,timestep = T],
            G[plant,timestep] <= gmax[plant]);

@constraint(Fax, EnergyBalance[node=NODES,timestep = T],
            sum(G[plant,timestep] for plant in plants_at_node[node])
            - demand[timestep, node] == INJ[node,timestep]);

@constraint(Fax, LineMax[line=LINES],
            sum(ptdf[line, node] * INJ[node,timestep] for node in NODES, timestep in T) <= l_max[line]
            );
@constraint(Fax, LineMin[line=LINES],
            sum(ptdf[line, node] * INJ[node,timestep] for node in NODES, timestep in T) >= -l_max[line]
            );

@constraint(Fax, sum(INJ[node,timestep] for node in NODES, timestep in T) == 0);

JuMP.optimize!(Fax)
JuMP.objective_value(Fax)
