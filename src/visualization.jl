# this file contains functions for visualizing the power grids
using Graphs
using CairoMakie
using GraphMakie
using CairoMakie.Colors
#### visualisation of square grid

### flow direction, definition of P_i:
### (plotted SimpleDiGraph: all horizontal lines point to the right, all vertical line point downwards)
# new: (compared to the old plotting flipped 180 degrees around the x-axis)
# if F_e < 0, flow rightwards or downwards (flow is in the same direction as arrow)
# if F_e > 0, flow leftwards or upwards (flow is in the opposite direction as arrow)
# P_i = 1 one unit of flow is generated
# P_i = -1 one unit of flow is consumed

### embed square grid
function set_vertex_locs(P) # sets vertex locations, example usage: gplot(g, locs_x, locs_y) while locs_x, locs_y = set_vertex_locs()
    N_vertices = length(P)
    N_side = sqrt(N_vertices)
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

### colour P_i
using Colors

# choice of colors: https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
function set_vertex_colors(P)
    nodefillc = RGB[] # initializing type
    for value in P
        if value == 1
            push!(nodefillc, colorant"grey85") # push! inserts items at end of collection # in TeX grey 95
        else
            push!(nodefillc, colorant"black") # in TeX grey 75
        end
    end
    nodefillc
end

function set_vertex_colors2(P)
    nodefillc = RGB[]
    for value in P
        if value == 1
            push!(nodefillc, colorant"grey95") # push! inserts items at end of collection # in TeX grey 95
        else
            push!(nodefillc, colorant"grey75") # in TeX grey 75
        end
    end
    nodefillc
end

function set_edge_colors(F::Array{Float64,1})
    edgefillc = RGB[]

    for value in F
        # push!(edgefillc, cgrad(:blues, Int(ceil(maximum(F)*100)), categorical = true)[Int(ceil((value*100)))])
        push!(edgefillc, cgrad(:blues, 101; categorical = true)[Int(ceil((value*100))+1)])

    end
    edgefillc
end



function set_gencon_colors(g,P)
    # loop over all edges
    edgefillc = RGB[]
    for i in 1:ne(g)
        if P[findfirst(isodd, B[:, i])] == P[findlast(isodd, B[:, i])]
            push!(edgefillc, colorant"grey85")

        else
            push!(edgefillc, colorant"black")
        end
    end
    edgefillc
end

function flows_colormap(g, P)
    F = flow(g, P)
    F_abs = abs.(F)
    locs_x, locs_y = set_vertex_locs(P)
    lay = X -> Point.(zip(locs_x,locs_y))

    # set edge colors
    edgecolors = set_edge_colors(F_abs)
    # set vertex colors
    vertex_colors = set_vertex_colors(P)

    fig = Figure() # creats Figure Object

    # determines position in figure
    ax1 = Axis(fig[1, 1], title ="Network", xgridvisible=false, ygridvisible=false)

    graphplot!(ax1, g, layout=lay, edge_width=F_abs*10,node_size=15, node_color=vertex_colors, edge_color=edgecolors)

    cbar = Colorbar(fig[1, 2], limits=(0,1), colormap = cgrad(:blues, 100),
    flipaxis = false, vertical = true,
    label = "flow units")
    cbar.ticks = 0:0.1:1

    fig
end


function compare_flows_colormap(g, P1, P2)
    F1 = flow(g, P1)
    F2 = flow(g, P2)
    F1_abs = abs.(F1)
    F2_abs = abs.(F2)

    locs_x, locs_y = set_vertex_locs(P1)
    lay = X -> Point.(zip(locs_x,locs_y))

    # set edge colors
    edgecolors1 = set_edge_colors(F1_abs)
    edgecolors2 = set_edge_colors(F1_abs)
    # set vertex colors
    vertex_colors1 = set_vertex_colors(P1)
    vertex_colors2 = set_vertex_colors(P2)


    fig = Figure() # creats Figure Object

    # determines position in figure
    ax1 = Axis(fig[1, 1], title ="Network 1", xgridvisible=false, ygridvisible=false)
    ax2 = Axis(fig[1, 2], title ="Network 2", xgridvisible=false, ygridvisible=false)

    # for manually adjusting postion and size of the suplots
    # see https://makie.juliaplots.org/stable/documentation/layoutables/
    # ax1 = Axis(fig, bbox = BBox(50, 375, 150, 475), title = "Network 1")
    # ax2 = Axis(fig, bbox = BBox(425, 750, 150, 475), title = "Network 2")

    # # adding edge and vertex labels
    # graphplot!(ax1, g, layout=lay, edge_width=F1_abs*10,node_size=15.0, node_color=vertex_colors1, edge_color=edgecolors1, elabels=string.(F1), elabels_textsize=12, nlabels=string.(P1), nlabels_textsize=12)
    # graphplot!(ax2, g, layout=lay, edge_width=F1_abs*10,node_size=15.0, node_color=vertex_colors2, edge_color=edgecolors2, elabels=string.(F2), elabels_textsize=12, nlabels=string.(P2), nlabels_textsize=12)
    graphplot!(ax1, g, layout=lay, edge_width=F1_abs*10,node_size=15.0, node_color=vertex_colors1, edge_color=edgecolors1)
    graphplot!(ax2, g, layout=lay, edge_width=F1_abs*10,node_size=15.0, node_color=vertex_colors2, edge_color=edgecolors2)

    # cbar = Colorbar(fig[2, 1:2], limits=(minimum(F),maximum(F)), colormap = cgrad(:blues, Int(ceil(maximum(F)*100))),
    cbar = Colorbar(fig[2, 1:2], limits=(0,1), colormap = cgrad(:blues, 100),
    flipaxis = false, vertical = false,
    label = "flow units")
    cbar.width = Relative(2/3)
    cbar.ticks = 0:0.1:1

    label_a = fig[1, 1, TopLeft()] = Label(fig, "A", textsize = 24,
    halign = :right)
    label_b = fig[1, 2, TopLeft()] = Label(fig, "B", textsize = 24,
    halign = :right)

    label_a.padding = (0, 6, 16, 0)
    label_b.padding = (0, 6, 16, 0)

    fig
end


function visualize_gencon(g, P)
    F = flow(g, P)
    F_abs = abs.(F)
    locs_x, locs_y = set_vertex_locs(P)
    lay = X -> Point.(zip(locs_x,locs_y))

    # set edge colors
    edgecolors = set_gencon_colors(g,P)
    # set vertex colors
    vertex_colors = set_vertex_colors(P)

    fig = Figure() # creats Figure Object

    # determines position in figure
    ax1 = Axis(fig[1, 1], title ="Network", xgridvisible=false, ygridvisible=false)

    # graphplot!(ax1, g, layout=lay, edge_width=F_abs*10,node_size=15.0, node_color=vertex_colors, edge_color=edgecolors, elabels=string.(F), elabels_textsize=12)
    graphplot!(ax1, g, layout=lay, edge_width=F_abs*10,node_size=15, node_color=vertex_colors, edge_color=edgecolors)

    fig
end




################################################################################
######################### DEPRECATED FUNCTIONS##################################
################################################################################

# # How to run deprecated functions using GraphPlot package:
#  - Do not load Graphs (do not execute `using Graphs`)
#  - load LightGraphs (using LightGraphs)
#  - load GraphPlot (using GraphPlot)
#  - load Colors (using Colors)
#
# old:
# ### flow direction, definition of P_i:
# ### (plotted SimpleDiGraph: all horizontal lines point to the right, all vertical line point downwards)
# # if F_e < 0, flow rightwards or downwards (flow is in the same direction as arrow)
# # if F_e > 0, flow leftwards or upwards (flow is in the opposite direction as arrow)
# # P_i = 1 one unit of flow is generated
# # P_i = -1 one unit of flow is consumed
#
# using GraphPlot
# ### show P_i in Graph
# # P = gen_rand_config(N_side)
# # nodelabel = P
# # gplot(g, nodelabel = P)
#
# ### show Flows F_e in Graph
# # F = flow(g, P)
# # g = gplot(g, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
#
# ### visualization of graph
# function visualize_graph(g, P)
#     F = flow(g, P)
#     locs_x, locs_y = set_vertex_locs(P)
#     nodefillc = set_vertex_colors(P)
#
#     gplot(g, locs_x, locs_y, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
#     #gplot(g, locs_x, locs_y, nodelabel = P, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
#     #gplot(g, locs_x, locs_y, nodefillc = nodefillc)
# end
#
# function visualize_graph_vlabel(g, P)
#     F = flow(g, P)
#     locs_x, locs_y = set_vertex_locs(P)
#     nodefillc = set_vertex_colors2(P)
#
#     #gplot(g, locs_x, locs_y, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
#     gplot(g, locs_x, locs_y, nodelabel = P, nodefillc = nodefillc, edgelabel = F, edgelabeldistx = 0, edgelabeldisty = 0)
#     #gplot(g, locs_x, locs_y, nodefillc = nodefillc)
# end
#
# function visualize_data(Data_loaded, P_rand_opt, Run_Nr) # for random grid: rand_opt=1 for optimized grid: rand_opt=4
#     P = Data_loaded[P_rand_opt][Run_Nr]
#     g = Data_loaded["Grid"]
#     visualize_graph(g, P)
# end
#
# ### visualize graph after line failure induced cascade
# # for evaluation of visualize_graph_after_linefailure_cascade: it must be line ⋜ m:
# #B = Array(incidence_matrix(g, oriented=true))
# #m = size(B)[2]
# function visualize_graph_after_linefailure_cascade(g, P, C, line)
#     g = linefailure!(g, line)
#     g = cascade!(g, P, C)
#     visualize_graph(g, P)
# end
#
# # save to pdf
# using Compose
# using Cairo
# draw(PDF("bla.pdf", 16cm, 16cm), gplot(g))
# # save to png
# draw(PNG("bla.png", 16cm, 16cm), gplot(g))
# # save to svg
# draw(SVG("bla.svg", 16cm, 16cm), gplot(g))
