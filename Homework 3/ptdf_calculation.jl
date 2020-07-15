using LinearAlgebra

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
    A = incedence
    Bd = Diagonal(b_vector)
    Bl = Bd*A # Line suceptance Matrix
    Bn = (A'*Bd)*A # Nodes suceptance Matrix

    B_inv = zeros(length(NODES),length(NODES))
    B_inv[setdiff(idx_nodes,idx_slack), setdiff(idx_nodes,idx_slack)] = inv(Bn[setdiff(idx_nodes,idx_slack), setdiff(idx_nodes,idx_slack)])
    PTDF = Bl*B_inv
    return PTDF
end

# Ehrenmann Smeers Example
function example_Ehrenmann()
    LINES = ["1-3","1-2","2-3","1-6","2-5","5-6","4-5","4-6"]
    NODES = ["N1", "N2", "N3", "N4", "N5", "N6"]
    SLACK = "N6"

    incedence = [1 0 -1 0 0 0; # 1 - 3
                 1 -1 0 0 0 0; # 1 - 2
                 0 1 -1 0 0 0; # 2 - 3
                 1 0 0 0 0 -1; # 1 - 6
                 0 1 0 0 -1 0; # 2 - 5
                 0 0 0 0 1 -1; # 5 - 6
                 0 0 0 1 -1 0; # 4 - 5-
                 0 0 0 1 0 -1; # 4 - 6
                 ]
    b_vector = [1;1;1;0.5;0.5;1;1;1]
    ptdf = calculate_ptdf(incedence, b_vector, NODES, SLACK)
    return ptdf
end
