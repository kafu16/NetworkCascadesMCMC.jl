#### building of square grid configuration
using LightGraphs, GraphPlot
using SyntheticNetworks


function gen_square_grid(N_side)
    N_vertices = N_side * N_side

    g = LightGraphs.SimpleGraph(N_vertices) # generates undirected graph

    # loop one: going through each row of square grid connecting each vertex with the vertex of its right
    for i in 1:N_side:N_vertices - N_side + 1
        for i in i:i + N_side - 2
            LightGraphs.add_edge!(g, i, i+1)
        end
    end

    # loop two: going through each column of square grid connecting each vertex with the vertex under it
    for i in 1:N_side
        for i in i:N_side:N_vertices - (2 * N_side) + i
            LightGraphs.add_edge!(g, i, i + N_side)
        end
    end
    g
end

function gen_periodic_square_grid(N_side)
    N_vertices = N_side * N_side

    g = LightGraphs.SimpleGraph(N_vertices) # generates undirected graph

    # loop one: going through each row of square grid connecting each vertex with the vertex of its right
    for i in 1:N_side:N_vertices - N_side + 1
        for i in i:i + N_side - 2
            LightGraphs.add_edge!(g, i, i+1)
        end
        LightGraphs.add_edge!(g, i, i + N_side - 1)
    end

    # loop two: going through each column of square grid connecting each vertex with the vertex under it
    for i in 1:N_side
        for i in i:N_side:N_vertices - (2 * N_side) + i
            LightGraphs.add_edge!(g, i, i + N_side)
        end
        LightGraphs.add_edge!(g, i, N_vertices - N_side + i)
    end
    g
end

""" generates power grid using SyntheticNetworks.jl and converts this generated
    grind into a SimpleGraph.
"""
function gen_rand_grid(rpg_seed, n, n0, p, q, r, s)
    Random.seed!(rpg_seed)
    u = 1 # parameter not implemented in SyntheticNetworks
    rpg = RandomPowerGrid(n, n0, p, q, r, s, u)
    generated_graph = generate_graph(rpg)
    A = LightGraphs.adjacency_matrix(generated_graph)
    g= LightGraphs.SimpleGraph(A)
end

# ### generation of non-embedded square grid (deprecated function)
# function gen_square_grid(N_side) # N_side: this number sqared gives number of vertices, for N_side > 2
#     N_vertices = N_side * N_side # N_vertices: number of vertices
#
#     g = LightGraphs.SimpleGraph(N_vertices) # generates undirected graph
#
#     # building of square grid
#     # loop one: going through each row of square grid (except last vertex of row, except last row) connecting
#     # each vertex with the vertex of its right and with the vertex under it
#     for i in 1:N_side:N_vertices - 2 * N_side +1
#         for i in i:i + N_side - 2
#             LightGraphs.add_edge!(g, i, i+1) # LightGraphs. because of "both Graphs and EmbeddedGraphs export "add_edge!"; uses of it in module Main must be qualified"
#             LightGraphs.add_edge!(g, i, i+N_side)
#         end
#     end
#
#     # loop two: going through last column of square grid connecting each vertex with the vertex under it
#     for i in N_side:N_side:N_vertices - N_side
#         LightGraphs.add_edge!(g, i, i + N_side)
#     end
#
#     # loop three: going through last row of square grid connecting each vertex with the vertex of its right
#     for i in N_vertices - N_side + 1:N_vertices -1
#         LightGraphs.add_edge!(g, i, i+1)
#     end
#     g
# end
