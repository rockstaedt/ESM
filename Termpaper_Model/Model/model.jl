###############################################################################
### Packages
###############################################################################

using JuMP
using Gurobi
using DataFrames
using LinearAlgebra
using CSV
using Plots
using Statistics

# Script to calculate PTDF Matrix
include("ptdf_calculation.jl")

###############################################################################
### Model Options
###############################################################################

# variable for including nuclear
nuclear_included = true

# variable for dc load flow
dc_load_flow = false

# time periods dictionary, value = path to dataframe
time_periods = Dict("1_winter_" => "data/48",
                    "2_spring_" => "",
                    "3_summer_" => "",
                    "4_autumn_" => ""
                    )

###############################################################################
### Data Import
###############################################################################

data_path = "data/48/"

L_df = joinpath(data_path,"lines_test.csv") |> CSV.read
N_df = joinpath(data_path,"nodes_test.csv") |> CSV.read
P_df = joinpath(data_path,"plants_test.csv") |> CSV.read

for (time_period, time_periode_path) in time_periods
    if time_periode_path != ""
        D_df = joinpath(time_periode_path,"demand_test.csv") |> CSV.read
        R_df = joinpath(time_periode_path,"avail_test.csv") |> CSV.read

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

        ## Power plants

        # dispatchable plants
        # check if nuclear energy is included
        if nuclear_included
            DISP = P_df[P_df.plant_type .== "conventional",:].index
        else
            DISP = P_df[(P_df.plant_type .== "conventional") .& (P_df.fuel .!= "uran"),:].index
        end
        # append hydro run on river
        append!(DISP, P_df[P_df.plant_type .== "hydro_ror",:].index)

        # nondispatchable plants
        # excluding "other_res" because no availability data
        NDISP = P_df[P_df.plant_type .== "wind onshore",:].index
        append!(NDISP, P_df[P_df.plant_type .== "wind offshore",:].index)
        append!(NDISP, P_df[P_df.plant_type .== "solar",:].index)

        PLANTS = vec(Array(select(P_df,"index")))

        # dataframe for DISP and NDISP

        # make empty copy
        DISP_PLANTS_df = copy(P_df)
        NDISP_PLANTS_df = copy(P_df)
        # delete information
        delete!(DISP_PLANTS_df,1:size(P_df,1))
        delete!(NDISP_PLANTS_df,1:size(P_df,1))
        # store information and reset index in Column1
        for (i, disp) in enumerate(DISP)
            f = P_df[P_df.index .== disp,:]
            # reset index -> be consistent with numeration in pandas
            f[!, "Column1"] .= i - 1
            append!(DISP_PLANTS_df,f)
        end
        for (i, ndisp) in enumerate(NDISP)
            f = P_df[P_df.index .== ndisp,:]
            # reset index -> be consistent with numeration in pandas
            f[!, "Column1"] .= i - 1
            append!(NDISP_PLANTS_df,f)
        end
        DISP_PLANTS_df
        NDISP_PLANTS_df

        ## Lines

        # dictionary to map line number with max flow
        line_number_to_maxflow = Dict(row.Column1 + 1 => row.maxflow for row in eachrow(L_df))

        # dictionary to map line index with max flow
        line_idx_to_maxflow = Dict(row.index => row.maxflow for row in eachrow(L_df))

        l_max = [line_number_to_maxflow[i] for i in 1:nrow(L_df)]

        ## Time

        T = T = 1:size(D_df, 1) |> collect

        ## Demand

        demand = D_df[:, 3:size(D_df,2)]

        ## Plants/node dictionary

        # dictionary to map disp plants to node
        disp_plants_at_node = Dict()
        for row in eachrow(DISP_PLANTS_df)
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
        for row in eachrow(NDISP_PLANTS_df)
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

        ## Marginal costs and max generation for conventional power plants

        # map marginal costs to conventional power plants
        mc_c = Dict(row.index => row.mc_el for row in eachrow(DISP_PLANTS_df))
        # map max of generation of conventional power plants
        gmax_c = Dict(row.index => row.g_max for row in eachrow(DISP_PLANTS_df))

        ##res feed dictionary for all timesteps and non dispatchable

        res_feed_in = Dict()
        for row in eachrow(NDISP_PLANTS_df)
                for t in T
                    res_feed_in[row.index,t] = row.g_max * R_df[t,row.index]
                end
        end
        res_feed_in

        ## Dictionaries for evaluation

        # map disp plant index to plant number
        disp_plant_idx_to_number = Dict(row.index => row.Column1+1 for row in eachrow(DISP_PLANTS_df))
        # dictionary to map disp plant number with node index
        disp_plant_number_to_idx = Dict(value => key for (key, value) in disp_plant_idx_to_number)

        # map ndisp plant index to plant number
        ndisp_plant_idx_to_number = Dict(row.index => row.Column1+1 for row in eachrow(NDISP_PLANTS_df))
        # dictionary to map disp plant number with node index
        ndisp_plant_number_to_idx = Dict(value => key for (key, value) in ndisp_plant_idx_to_number)

        # map plant index with fuel type
        plant_idx_to_fuel = Dict()
        for row in eachrow(P_df)
            plant_idx_to_fuel[row.index] = row.fuel
        end
        plant_idx_to_fuel

        ###############################################################################
        ### Model DISPATCH
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

        @constraint(m, Slack,
                    sum(INJ[node,t] for node in NODES, t in T) == 0
                    );

        if !dc_load_flow
            JuMP.optimize!(m)
        end

        ###############################################################################
        ### Model DC LOAD FLOW
        ###############################################################################

        if dc_load_flow
            @constraint(m, LineMax[line=LINES, t=T],
                        sum(ptdf[line_idx_to_number[line],node_idx_to_number[node]] * INJ[node,t] for node in NODES)
                        <= line_idx_to_maxflow[line]
                        );

            @constraint(m, LineMin[line=LINES, t=T],
                        sum(ptdf[line_idx_to_number[line],node_idx_to_number[node]] * INJ[node,t] for node in NODES)
                        >= -line_idx_to_maxflow[line]
                        );

            JuMP.optimize!(m)
        end

        ###############################################################################
        ### Results
        ###############################################################################

        # check for optimal solution
        if termination_status(m) == MOI.OPTIMAL
            println()
            println("##############################################################")
            println()
            println("Optimal solution found.")
            println()

            # check for model type and nuclear and set path and model_type
            if dc_load_flow & nuclear_included
                model_type = "DC Load Flow with nuclear"
                export_path = "export_files/dc_load_flow_nuclear/"
            elseif dc_load_flow & !nuclear_included
                model_type = "DC Load Flow without nuclear"
                export_path = "export_files/dc_load_flow/"
            elseif !dc_load_flow & nuclear_included
                model_type = "Dispatch with nuclear"
                export_path = "export_files/dispatch_nuclear/"
            elseif !dc_load_flow & !nuclear_included
                model_type = "Dispatch without nuclear"
                export_path = "export_files/dispatch/"
            end

            # create path
            mkpath(export_path)

            println("--------------------------------------------------------------")
            println("Results for model: ", model_type, " - Period: ", time_period)
            println("--------------------------------------------------------------")
            println("--------------------------------------------------------------")
            println()

            # Generation costs
            #---------------------------------------------------------------------------
            println("Cost of generation: ", JuMP.objective_value(m))
            println()

            # Generation dataframe
            #---------------------------------------------------------------------------

            println("--> Generating generation dataframe...")
            println()

            disp_generation = value.(G).data
            ndisp_generation = [res_feed_in[ndisp, t] for ndisp in NDISP, t in T]

            # build generation dataframe by fuel type
            generation_df_fuel = DataFrame(timestep = T)
            fuel_types = unique(P_df, "fuel").fuel
            for fuel_type in fuel_types
                generation_df_fuel[!, fuel_type] .= 0.0
            end
            # dispatchable generation
            # rows: plant_number, columns: timesteps
            for timestep in 1:size(disp_generation, 2)
                for plant_number in 1:size(disp_generation, 1)
                    generation_df_fuel[timestep, plant_idx_to_fuel[disp_plant_number_to_idx[plant_number]]] += disp_generation[plant_number, timestep]
                end
            end
            # dispatchable generation
            # rows: plant_number, columns: timesteps
            for timestep in 1:size(ndisp_generation, 2)
                for plant_number in 1:size(ndisp_generation, 1)
                    generation_df_fuel[timestep, plant_idx_to_fuel[ndisp_plant_number_to_idx[plant_number]]] += ndisp_generation[plant_number, timestep]
                end
            end

            CSV.write(joinpath(export_path, string(time_period, "generation_by_fuel.csv")), generation_df_fuel)
            println("--> Exporting of generation file...")
            println()

            # Curtailment
            #---------------------------------------------------------------------------

            println("--> Calculating curtailment...")
            println()

            curtailment = value.(CU).data

            for i in 1:size(curtailment, 1)
                for j in 1:size(curtailment, 2)
                    if curtailment[i,j] > 0
                        println("Curtailment in t: ", j, " and node: ", node_number_to_idx[i])
                        println("Amount: ", curtailment[i,j])
                    end
                end
            end

            println()

            # Consumption
            #---------------------------------------------------------------------------

            println("--> Generating consumption dataframe...")
            println()

            consumption_df = DataFrame(timestep = T)
            consumption_df[!, "demand"] .= 0.0

            # add demand
            demand_matrix = convert(Matrix, demand)
            for row in 1:size(demand_matrix, 1)
                consumption_df[row, "demand"] = -sum(demand_matrix[row,:])
            end
            # add curtailment
            consumption_df[!, "curtailment"] .= 0.0
            curtailment_matrix = curtailment'
            for row in 1:size(curtailment_matrix, 1)
                consumption_df[row, "curtailment"] = -sum(curtailment_matrix[row,:])
            end

            CSV.write(joinpath(export_path, string(time_period, "consumption.csv")), consumption_df)
            println("--> Exporting consumption file...")
            println()

            # Injections
            #---------------------------------------------------------------------------

            injections = value.(INJ).data

            # Line results
            #---------------------------------------------------------------------------

            line_results_t = Dict()

            for t in T
                df = DataFrame(Lines=LINES,
                                Flow=ptdf*value.(INJ).data[:,t],
                                Line_Capacity = [line_idx_to_maxflow[line] for line in LINES]
                                )
                df[!, "Capacity"] .= abs.(df.Flow ./ df.Line_Capacity)
                line_results_t[t] = df
            end

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

            # Price
            #---------------------------------------------------------------------------

            println("--> Calculating price...")
            println()

            # unit: €/MWh
            price=dual.(EnergyBalance).data

            global min_price = Dict("price" => price[1,1], "node" => 1, "timestep" => 1)
            global max_price = Dict("price" => price[1,1], "node" => 1, "timestep" => 1)
            for i in 1:size(price, 1)
                global min_price
                global max_price
                for j in 1:size(price, 2)
                    if price[i,j] < min_price["price"]
                        min_price["price"] = price[i,j]
                        min_price["node"] = i
                        min_price["timestep"] = j
                    end
                    if price[i,j] > max_price["price"]
                        max_price["price"] = price[i,j]
                        max_price["node"] = i
                        max_price["timestep"] = j
                    end
                end
            end

            println("Minium price = ", min_price["price"], " in timestep: ", min_price["timestep"], " and node: ", node_number_to_idx[min_price["node"]])
            println()

            println("Maximum price = ", max_price["price"], " in timestep: ", max_price["timestep"], " and node: ", node_number_to_idx[max_price["node"]])
            println()

            println("Average price = ", mean(price))
            println()

            println("##############################################################")
            println()
        else
            println()
            println("##############################################################")
            println()
            println("Optimal solution not found. Please check!")
            println()
            println("##############################################################")
            println()
        end
    end
end
## Snippets

# generation_t2 = DataFrame(Plants = PLANTS,
#                         t = "t2",
#                         DISP = value.(G).data[:,2]
#                         )

#generation = vcat(generation_t1, generation_t2)


#CSV.write(joinpath(export_path, "generation.csv"), generation)
