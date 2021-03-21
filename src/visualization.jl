# this file contains functions for visualizing the power grids

#### visualisation of square grid

### embed square grid
function set_vertex_locs(N_side) # sets vertex locations, example usage: gplot(g, locs_x, locs_y) while locs_x, locs_y = set_vertex_locs()
    N_vertices = N_side * N_side
    locs_x = zeros(Int64, N_vertices) # generation of arrays
    locs_y = zeros(Int64, N_vertices)

    # generation of locs_x
    k = 1
    j = 1
    for i in 1:N_side
        for j in 1:N_side
            locs_x[k] = j
            k = k + 1
            j = j + 1
        end
        j = 1
    end

    # generation of locs_y
    k = 1
    for i in 1:N_side
        for j in 1:N_side
            locs_y[k] = i
            k = k + 1
        end
    end

    locs_x = convert(Vector{Float64}, locs_x) # conversion from Vector{Int64} to Vector{Float64}
    locs_y = convert(Vector{Float64}, locs_y)
    locs_x, locs_y
end

### show P_i in Graph
using Graphs

# P = gen_rand_config(N_side)
# nodelabel = P
# gplot(g, nodelabel = P)

### colour P_i
using Colors

# choice of colours: https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
function set_colours(P)
    nodefillc = [] # generates empty vector of type 'Any'
    for value in P
        if value == 1
            push!(nodefillc, colorant"grey85") # push! inserts items at end of collection # in TeX grey 95
        else
            push!(nodefillc, colorant"black") # in TeX grey 75
        end
    end
    nodefillc
end

function set_colours2(P)
    nodefillc = [] # generates empty vector of type 'Any'
    for value in P
        if value == 1
            push!(nodefillc, colorant"grey95") # push! inserts items at end of collection # in TeX grey 95
        else
            push!(nodefillc, colorant"grey75") # in TeX grey 75
        end
    end
    nodefillc
end
### show Flows F_e in Graph
# F = flow(g, P)
# g = gplot(g, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)

### visualization of graph
function visualize_graph(g, P, N_side)
    F = flow(g, P)
    locs_x, locs_y = set_vertex_locs(N_side)
    nodefillc = set_colours(P)

    gplot(g, locs_x, locs_y, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
    #gplot(g, locs_x, locs_y, nodelabel = P, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
    #gplot(g, locs_x, locs_y, nodefillc = nodefillc)
end

function visualize_graph_vlabel(g, P, N_side)
    F = flow(g, P)
    locs_x, locs_y = set_vertex_locs(N_side)
    nodefillc = set_colours2(P)

    #gplot(g, locs_x, locs_y, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
    gplot(g, locs_x, locs_y, nodelabel = P, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
    #gplot(g, locs_x, locs_y, nodefillc = nodefillc)
end

### flow direction, definition of P_i:
### (plotted SimpleDiGraph: all horizontal lines point to the right, all vertical line point downwards)
# if F_e < 0, flow rightwards or downwards (flow is in the same direction as arrow)
# if F_e > 0, flow leftwards or upwards (flow is in the opposite direction as arrow)
# P_i = 1 one unit of flow is generated
# P_i = -1 one unit of flow is consumed

### visualize graph after line failure induced cascade
# for evaluation of visualize_graph_after_linefailure_cascade: it must be line ⋜ m:
#B = Array(incidence_matrix(g, oriented=true))
#m = size(B)[2]
function visualize_graph_after_linefailure_cascade(P, C, N_side, line)
    g = gen_square_grid(N_side)
    g = linefailure!(g, line)
    g = cascade!(g, P, C)
    visualize_graph(g, P, N_side)
end
