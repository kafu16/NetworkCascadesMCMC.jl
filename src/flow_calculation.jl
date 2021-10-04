### flow calculation
using IterativeSolvers
using LightGraphs


# N_vertices: number of vertices
# m: number of edges
function flow(g, P) # calculates the flow of graph g for a certain configuration P
    B = incidence_matrix(g, Float64, oriented=true) # nxm oriented incidence matrix B
    # node-edge matrix, indexed by [v, i], i is in 1:ne(g),
    # indexing edge e. Directed graphs: -1 indicates src(e) == v, 1 indicates dst(e) == v.
    # indexing of the edges e does not correspond to sequence of building edges above!!!
    # indexing of edges: https://github.com/JuliaGraphs/LightGraphs.jl/blob/216f5ffa77860d4c39b8d05fe2197d0af75a4241/src/linalg/spectral.jl#L129-L142

    F = lsqr(B, P) # solving BF = P for F by returning minimal norm solution: https://juliamath.github.io/IterativeSolvers.jl/dev/linear_systems/lsqr/
    F_rounded = round.(F; digits = 2)
end

# function flow(g, P) # calculates the flow of graph g for a certain configuration P
#     L = Array(LightGraphs.laplacian_matrix(g)) # N_vertices x N_vertices Laplacian Matrix L
#     # Array(L) converts sparse matrix to dense matrix
#     # https://docs.julialang.org/en/v1/stdlib/SparseArrays/index.html#Correspondence-of-dense-and-sparse-methods-1
#
#     B = Array(incidence_matrix(g, oriented=true))
#     R = pinv(L)
#
#     φ = R * P # φ_i: potential of node i
#     F = transpose(B) * φ # F_e: Flow of edge e, F: Vector of all flows
#     F_rounded = round.(F; digits = 2)
# end


# function flow(g, P) # calculates the flow of graph g for a certain configuration P
#     B = Array(incidence_matrix(g, oriented=true)) # nxm oriented incidence matrix B
#     # node-edge matrix, indexed by [v, i], i is in 1:ne(g),
#     # indexing edge e. Directed graphs: -1 indicates src(e) == v, 1 indicates dst(e) == v.
#     # indexing of the edges e does not correspond to sequence of building edges above!!!
#     # indexing of edges: https://github.com/JuliaGraphs/LightGraphs.jl/blob/216f5ffa77860d4c39b8d05fe2197d0af75a4241/src/linalg/spectral.jl#L129-L142
#     R = pinv(B)
#
#     F = R * P
#     F_rounded = round.(F; digits = 2)
# end
