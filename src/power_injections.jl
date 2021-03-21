# net power injections P_i

#### Building of an array whose entries represent nodes of square grid.
# Vertex one is entry one of array counting until the end of first line
# of square grid, continuing to count in second line and so on.
# Above built graph and array are (state: 25.09.2019) non-interdependent.

#  initial configuration: flow network with random placements nodes having inflow = 1
# or outflow = -1 representing consumers and generators

using Random

function gen_rand_config(N_side) # generates random configuration P (see below)
    N_vertices = N_side * N_side
    P = ones(Float64, N_vertices) # P_i: net inflow or outflow at vertex i
    P[1:2:end] .= P[1:2:end] * -1. # every second element is changed from 1 to -1
    #### ToDo evtl. '.' entfernen, weil nicht einheitlich
    shuffle!(P) # randomly permutes entries of array
end

function gen_stable_config(g, N_side, C) # generation of initial stable random configuration P
#### ToDo build: return error when stable config is not possible due to a too low C
    P = gen_rand_config(N_side)
    F = flow(g, P)
    while maximum(abs.(F)) > C
        P = gen_rand_config(N_side)
        F = flow(g, P)
    end
    P
end
