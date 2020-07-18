using JuMP
using Ipopt ## Needed for Quadratic Programs
using DataFrames

include("ptdf_calculation.jl")

NODES = ["N1", "N2", "N3", "N4", "N5", "N6"]
SLACK = "N6"
DEM_NODES = ["N3","N5","N6"]
SUP_NODES = ["N1","N2","N4"]


LINES = ["1-3","1-2","2-3","1-6","2-5","5-6","4-5","4-6"]

## Topology for Ehrenmann Smeers Network

incedence = [1 0 -1 0 0 0; # 1 - 3
             1 -1 0 0 0 0; # 1 - 2
             0 1 -1 0 0 0; # 2 - 3
             1 0 0 0 0 -1; # 1 - 6
             0 1 0 0 -1 0; # 2 - 5
             0 0 0 0 1 -1; # 5 - 6
             0 0 0 1 -1 0; # 4 - 5-
             0 0 0 1 0 -1; # 4 - 6
            ]
## Corresponding susceptance values
# Paper states impedance to be 1, except for two lines, where they are 2
b_vector = [1;1;1;0.5;0.5;1;1;1]

# Calculate PTDF Matrix from function
PTDF = calculate_ptdf(incedence, b_vector, NODES, SLACK)

# Put into Dictionary, so that we can access it via ptdf[line, node]
ptdf = Dict((l,n) => PTDF[l_idx, n_idx] for (n_idx, n) in enumerate(NODES),
                                            (l_idx,l) in enumerate(LINES))

a = Dict("N1" => 10.0,
         "N2" => 15.0,
         "N3" => 37.5,
         "N4" => 42.5,
         "N5" => 75.0,
         "N6" => 80.0)

b = Dict("N1" => 0.05,
         "N2" => 0.05,
         "N3" => -0.05,
         "N4" => 0.025,
         "N5" => -0.1,
         "N6" => -0.1)

# from paper Chao and Peck (1998)
lmax_dict = Dict("1-6" => 200, "2-5" => 250)

Ehrenmann = Model(Ipopt.Optimizer)
@variables Ehrenmann begin
    Q[NODES] >= 0
    INJ[NODES]
end

# maximation of total welfare (demand function - supply function or consumer surplus + producer surplus)
@objective(Ehrenmann, Max,
    sum((a[dem_node] + 0.5*b[dem_node]*Q[dem_node])*Q[dem_node] for dem_node in DEM_NODES)
    - sum((a[sup_node] + 0.5*b[sup_node]*Q[sup_node])*Q[sup_node] for sup_node in SUP_NODES)
    );

@constraint(Ehrenmann, EnergyBalance[node=NODES],
    sum(Q[sup_node] for sup_node in intersect([node], SUP_NODES)) -
    sum(Q[dem_node] for dem_node in intersect([node], DEM_NODES)) == INJ[node]
    );

@constraint(Ehrenmann, Slack,
    sum(INJ[node] for node in NODES) == 0
    );


## Implement LineMax and LineMin for lines "1-6" and "2-5"
# Task 2 b)

# Implement Zonal Price Constraint:
## price_n for n in Price Zone == Price of that Price Zone
# Task 2 c)

JuMP.optimize!(Ehrenmann)

println("")
println("Objective Value: $(JuMP.objective_value(Ehrenmann))")

nodal_results = DataFrame(Nodes=NODES,
                          Quantity=value.(Q).data,
                          Injections=value.(INJ).data,
                          Price_Dual=dual.(EnergyBalance).data)

nodal_results[!, :Price_inv] = [a[node] + b[node]*sum(nodal_results[nodal_results[:, :Nodes] .== node, :Quantity]) for node in NODES]

line_cap = [(line in keys(lmax_dict) ? lmax_dict[line] : 0) for line in LINES]

line_results = DataFrame(Lines=LINES,
                        Flow=PTDF*value.(INJ).data,
                        Line_Capacity = line_cap)

println(nodal_results)
println(line_results)
