# net power injections P_i

using Random

################################ square grids ##################################
# Building of an array whose entries represent nodes of square grid.
# Vertex one is entry one of array counting until the end of first line
# of square grid, continuing to count in second line and so on.

#  initial configuration: flow network with random placements nodes having inflow = 1
# or outflow = -1 representing consumers and generators

function gen_rand_config(N_vertices) # generates random configuration P (see below)
    P = ones(Float64, N_vertices) # P_i: net inflow or outflow at vertex i
    P[1:2:end] .= P[1:2:end] * -1. # every second element is changed from 1 to -1
    shuffle!(P) # randomly permutes entries of array
end

function gen_stable_config(g, C) # generation of initial stable random configuration P
#### ToDo build: return error when stable config is not possible due to a too low C
    N_vertices = nv(g)
    P = gen_rand_config(N_vertices)
    F = flow(g, P)
    while maximum(abs.(F)) > C
        P = gen_rand_config(N_vertices)
        F = flow(g, P)
    end
    P
end

function gen_multiple_stable_configs(g, C, N_runs)
    P_inits = []
    for i in 1:N_runs
        P_init = gen_stable_config(g, C)
        push!(P_inits, P_init)
    end
    P_inits
end
