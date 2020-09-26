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


###############################################################################
### Evaluation
###############################################################################


function evaluate_model(model, value, model_type)

    # check for optimal solution
    if termination_status(m) == MOI.OPTIMAL
        println("##############################################################")
        println()
        println("Optimal solution found.")
        println()

        # check for model type
        if model_type == "dispatch"
            export_path = "export_files/disptach/"
            filename_extension = "-dispatch"
        else
            export_path = "export_files/dc_load_flow//"
            filename_extension = "-dcloadflow"
        end

        println("--------------------------------------------------------------")
        println("Results for model: ", model_type)
        println("--------------------------------------------------------------")
        println()

        ########################################################################
        # Generation costs
        ########################################################################
        println("Cost of generation: ", JuMP.objective_value(m))

        ########################################################################
        # Generation dataframe
        ########################################################################

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

        generation_df_fuel

        CSV.write(joinpath(export_path, "generation_by_fuel.csv"), generation_df_fuel)



        println("##############################################################")
    else
        println("##############################################################")
        println()
        println("Optimal solution not found. Please check!")
        println()
        println("##############################################################")
    end

end
