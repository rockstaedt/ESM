###############################################################################
### Packages
###############################################################################

using JuMP
using Gurobi
using DataFrames
using LinearAlgebra
using CSV
using Plots


# Script to calculate PTDF Matrix
include("ptdf_calculation.jl")

###############################################################################
### Data Import
###############################################################################

data_path = "data"

L_df = joinpath(data_path,"lines_test.csv") |> CSV.read
N_df = joinpath(data_path,"nodes_test.csv") |> CSV.read
P_df = joinpath(data_path,"plants_test.csv") |> CSV.read
D_df = joinpath(data_path,"demand_test.csv") |> CSV.read
R_df = joinpath(data_path,"avail_test.csv") |> CSV.read

###############################################################################
### Data Preprocessing
###############################################################################

## Create PTDF Matrix

# empty incidence matrix
A = zeros(nrow(L_df),nrow(N_df))

# dictionary to map node index with a number for a fixed order
node_idx_to_number = Dict(row.index => row.Column1+1 for row in eachrow(N_df))

# dictionary to map node number with node index
node_number_to_idx = Dict(value => key for (key, value) in node_idx_to_number)

# list of all nodes
NODES = [node_number_to_idx[i] for i in 1:nrow(N_df)]

# dictionary to map line index with a number for a fixed order
line_idx_to_number = Dict(row.index => row.Column1+1 for row in eachrow(L_df))

# dictionary to map node number with node index
line_number_to_idx = Dict(value => key for (key, value) in line_idx_to_number)

# list of all lines
LINES = [line_number_to_idx[i] for i in 1:nrow(L_df)]

# dictionary to map line number with susceptance
line_number_to_sus = Dict(row.Column1 + 1 => row.b for row in eachrow(L_df))

# susceptance vector
b = [line_number_to_sus[i] for i in 1:nrow(L_df)]

# slack node
SLACK = NODES[1]

# calculation of incidence matrix A
for row in eachrow(L_df)
    line_number = row.Column1+1
    node_start = row.node_i
    node_end = row.node_j
    A[line_number,node_idx_to_number[node_start]] = 1
    A[line_number,node_idx_to_number[node_end]] = -1
end
A

ptdf = calculate_ptdf(A,b,NODES,SLACK)

# dispatchable plants
DISP = P_df[P_df.plant_type .== "conventional",:].index
append!(DISP, P_df[P_df.plant_type .== "hydro_ror",:].index)
# nondispatchable plants
# excluding "other_res" because no availability data
NDISP = P_df[P_df.plant_type .== "wind onshore",:].index
append!(NDISP, P_df[P_df.plant_type .== "wind offshore",:].index)
append!(NDISP, P_df[P_df.plant_type .== "solar",:].index)

PLANTS = vec(Array(select(P_df,"index")))

# dictionary to map line number with max flow
line_number_to_maxflow = Dict(row.Column1 + 1 => row.maxflow for row in eachrow(L_df))

# dictionary to map line index with max flow
line_idx_to_maxflow = Dict(row.index => row.maxflow for row in eachrow(L_df))

l_max = [line_number_to_maxflow[i] for i in 1:nrow(L_df)]


T = T = 1:size(D_df, 1) |> collect
demand = D_df[:, 3:size(D_df,2)] #demand dataframe


#df for disp plants

#make empty copy
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
NONDISP_PLANTS

# dictionary to map disp plants to node
disp_plants_at_node = Dict()
for row in eachrow(DISP_PLANTS)
    if !haskey(disp_plants_at_node,row.node)
        disp_plants_at_node[row.node]=[row.index]
    else
        push!(disp_plants_at_node[row.node],row.index)
    end
end
# check if dictionary for disp plants contains all nodes
# -> it can happen that a node only contains renewables or conventional powerplants
if length(disp_plants_at_node) < length(NODES)
    for node in NODES
        if !haskey(disp_plants_at_node, node)
            disp_plants_at_node[node] = []
        end
    end
end
disp_plants_at_node

# dictionary to map ndisp plants to node
ndisp_plants_at_node = Dict()
for row in eachrow(NONDISP_PLANTS)
    if !haskey(ndisp_plants_at_node,row.node)
        ndisp_plants_at_node[row.node]=[row.index]
    else
        push!(ndisp_plants_at_node[row.node],row.index)
    end
end
# check if dictionary for non disp plants contains all nodes
# -> it can happen that a node only contains renewables or conventional powerplants
if length(ndisp_plants_at_node) < length(NODES)
    for node in NODES
        if !haskey(ndisp_plants_at_node, node)
            ndisp_plants_at_node[node] = []
        end
    end
end
ndisp_plants_at_node

# map marginal costs to conventional power plants
mc_c = Dict(row.index => row.mc_el for row in eachrow(DISP_PLANTS))
# map max of generation of conventional power plants
gmax_c = Dict(row.index => row.g_max for row in eachrow(DISP_PLANTS))

#plants_at_node[SLACK]=[]
#to enter slack
#Kopie = copy(P_df)
#push!(Kopie,["0" 14 "0" 0 0 "0" 0])

#gmax_r für alle nodes, ts
res_feed_in = Dict()
for row in eachrow(NONDISP_PLANTS)
        for t in T
            res_feed_in[row.index,t] = row.g_max * R_df[t,row.index]
        end
end
res_feed_in

###############################################################################
### Model
###############################################################################

m = Model(Gurobi.Optimizer)

@variables m begin
    G[DISP,T] >= 0
    INJ[NODES,T]
    CU[NODES, T] >= 0
end

@objective(m, Min, sum(mc_c[disp] * G[disp,t] for disp in DISP, t in T));

@constraint(m, Max_Generation[disp=DISP,t=T],
            G[disp,t] <= gmax_c[disp]);

@constraint(m, EnergyBalance[node=NODES,t=T],
            sum(G[plant,t] for plant in disp_plants_at_node[node])
            + sum(res_feed_in[ndisp, t] for ndisp in ndisp_plants_at_node[node])
            - demand[t,node]
            - CU[node, t]
            == INJ[node,t]);

@constraint(m, LineMax[line=LINES, t=T],
            sum(ptdf[line_idx_to_number[line],node_idx_to_number[node]] * INJ[node,t] for node in NODES)
            <= line_idx_to_maxflow[line]
            );

@constraint(m, LineMin[line=LINES, t=T],
            sum(ptdf[line_idx_to_number[line],node_idx_to_number[node]] * INJ[node,t] for node in NODES)
            >= -line_idx_to_maxflow[line]
            );

@constraint(m, Slack,
            sum(INJ[node,t] for node in NODES, t in T) == 0
            );

JuMP.optimize!(m)
JuMP.objective_value(m)


###############################################################################
### Results
###############################################################################

export_path = "export_files"

generation = value.(G).data
injections = value.(INJ).data
curtailment = value.(CU).data

for i in 1:size(curtailment, 1)
    for j in 1:size(curtailment, 2)
        if curtailment[i,j] > 0
            println("Curtailment in t: ", j, " and node: ", node_number_to_idx[i])
            println("Amount: ", curtailment[i,j])
        end
    end
end

line_results_t = Dict()

for t in T
    df = DataFrame(Lines=LINES,
                    Flow=ptdf*value.(INJ).data[:,t],
                    Line_Capacity = [line_idx_to_maxflow[line] for line in LINES]
                    )
    df[!, "Capacity"] .= abs.(df.Flow ./ df.Line_Capacity)
    line_results_t[t] = df
end

line_results_t

# for (key, value) in line_results_t
#     for row in eachrow(value)
#         if row.Capacity >= 0.99
#             print("Capacity reaches nearly limit on line: ")
#             println(row.Lines)
#             print("At timeslot: ")
#             println(key)
#             print("Capacity: ")
#             println(row.Capacity)
#             println("-----")
#         end
#     end
# end

# generation_t2 = DataFrame(Plants = PLANTS,
#                         t = "t2",
#                         DISP = value.(G).data[:,2]
#                         )

#generation = vcat(generation_t1, generation_t2)


#CSV.write(joinpath(export_path, "generation.csv"), generation)

# unit:€/MWh
price=dual.(EnergyBalance).data

min_price = price[1,1]
max_price = price[1,1]
for i in 1:size(price, 1)
    global min_price
    global max_price
    for j in 1:size(price, 2)
        if price[i,j] < min_price
            min_price = price[i,j]
        end
        if price[i,j] > max_price
            max_price = price[i,j]
        end
    end
end
min_price
max_price

using Statistics

mean(price)
