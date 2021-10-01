################################################################################
# grid operations (funtions that change/refer to  grid topology) used for SA ###
################################################################################

using LightGraphs

### failure/removal of edge that causes cascade
function linefailure!(g::LightGraphs.AbstractGraph, i::Integer) # edge i is removed
    B = Array(incidence_matrix(g, oriented=true))
    rem_edge!(g, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
    g
end


### cascade
function cascade!(g::LightGraphs.AbstractGraph, P::Array{Float64,1}, C::AbstractFloat) # removes all edges whose flow is bigger than C -> calculates flow again and so on until all flows are smaller than C
    F = flow(g, P)
    while maximum(abs.(F)) > C

        B = Array(incidence_matrix(g, oriented=true))
        m = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
        # it is ne(g) = size(B)[2], where ne() give number of edges
        for i in 1:m # in this loop all edges that carry flow bigger than certain value/ threshold are being removed
            if abs(F[i]) > C
                rem_edge!(g, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
                # B[:, i] gives elements of i-th column of B, isodd accesses only the odd elements of these element
                # findfirst()/ findlast() access first/ last odd element.
                # incidence matrix B comprises only zeros (even elements) and one 1 and one -1 for the two vertices that are connected
                # to respective edge
            end
        end
        F = flow(g, P)
    end
    g
end


### measurements of energy() (see below)
function biggest_component(g::LightGraphs.AbstractGraph) # gives biggest connected component of graph g
    maximum(length.(LightGraphs.connected_components(g))) # connected_components() returns vector whose components contain connected vertices
    # so the length of a component is the size of a connected graph
    # so here the biggest connected graph out of all vertices is chosen
end


################################################################################
############################# energy function(s) ###############################
################################################################################
using Statistics

function energy(g_init::LightGraphs.AbstractGraph, P::Array{Float64,1}, C::AbstractFloat) # calculates energy of step k, C: threshold that marks line failure,
    g = copy(g_init)
    B = Array(incidence_matrix(g, oriented=true))
    m = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
    linefailure_indizes = collect(1:m) # collect() collects all values in the range 1:m in an array, here all edges are numbered in an array
    linefailure_indizes = shuffle!(linefailure_indizes) # positions of edges in linefailure_indizes is randomly permuted by shuffle!()
    # randomness is redundant unless m in following for-loop is replaced by a number smaller than m

    # # N_removals (optional argument): number of edge removals for approximation
    # # (not implemented)
    # if N_removals > 0
    #     N = N_removals
    # else
    #     N = m
    # end
    # # Für Näherung, d.h. man zieht nur N edges anstatt alle, muss man das nachfolgenden m durch N ersetzen

    G = [ ]
    for i in 1:m # for loop for randomly chosen N_vertices linefailures
        g = linefailure!(g, linefailure_indizes[i])
        g = cascade!(g, P, C)

        G = append!(G, biggest_component(g))
        g = copy(g_init)
    end

    G_av = mean(G)
    G, G_av # this way two values in a tuple are returned by a function
end


################################################################################
########################## Monte Carlo step functions ##########################
################################################################################

### swap of configurations
function swap!(P::Array{Float64,1}) # swaps randomly chosen generator-consumer pair
    N_vertices = length(P)
    A = rand(1:N_vertices,1)
    B = rand(1:N_vertices,1)
    while A == B || P[A] == P[B]
        A = rand(1:N_vertices,1)
        B = rand(1:N_vertices,1)
    end
    P[A], P[B] = P[B], P[A] # swaps neighbors
    P
end


# generation of stable swapped configuration
function stable_swapped_config!(g::LightGraphs.AbstractGraph, P::Array{Float64,1}, C::AbstractFloat)
# to avoid calculation steps, the input configuration should be stable as the stable_swapped_config() only permutes
# one generator-consumer pair. so having an unstable configuration as input will probably take more steps than
# first generating a stable configuration by gen_stable__square_config() and then apply stable_swapped_config()
#### ToDo build: return error when stable config is not possible due to a too low C
    P_stable_old = copy(P)
    # by using P_stable_old it is made sure that a the swapped configuration differs only by one permutation
    # otherwise by the following while-loop in each iteration of the loop the newly generated configuration would be swapped
    # again to obtain a stable configuration. This can not be done by P_stable_old = P because of variable assignment. In the latter case
    # P_stable_old would equal P all the time, even when the value of P is changed
    P = swap!(P)
    F = flow(g, P)
    #### ToDo: implement without `global`.
    # max_iterations = 10000 # hardcoded number of iteration after that loop is excited
    # global i = 1
    while maximum(abs.(F)) > C
        P = swap!(P_stable_old)
        F = flow(g, P)
        #### ToDo
        # if i >= max_iterations
        #     error("ERROR: Maximum number of iterations for finding stable configuration reached.")
        # end
        # global i += 1
    end
    P
end


################################################################################
############### implementation of Simulated Annealing (SA) #####################
################################################################################

### key code of SA

# core function for simulated annealing
function sim_anneal(g::LightGraphs.AbstractGraph, P_init::Array{Float64,1}, C::AbstractFloat, annealing_schedule::Function, steps_per_temp::Integer, k_max::Integer) # k_max: setting number of computation steps
    # given an initial configuration P sim_anneal() tendentially finds a more stable configuration

    en = [ ]
    energy_init = energy(g, P_init, C)[2] # by [2] only the second value of tuple is returned (G_av)
    push!(en, energy_init)
    P = copy(P_init)

    for k in 1:k_max-1
        T = annealing_schedule(k,steps_per_temp) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
        P_old = copy(P) # for calculating the energy of "old" configuration
        P = stable_swapped_config!(g, P, C)
        energy_old = en[k]
        energy_new = energy(g, P, C)[2] # by [2] only the second value of tuple is returned (G_av)
        ΔE = energy_new - energy_old
        #### performance: let energy() calculate G_av only? Nope, G is only saving
        # to an array and is helpful for understanding the algorithm.

        if ΔE <= 0 # man könnte conditional auch umdrehen: if (ΔE <= 0 AND probability(ΔE, T) < rand())
                                                                 # P = P_old
            push!(en, energy_new)
        elseif probability(ΔE, T) > rand() # rand() gives random number element of [0,1]
            push!(en, energy_new)
        else
            P = P_old
            push!(en, energy_old)
        end

    end
    P, en
end


function multiple_sim_anneal(filepath::String, g::LightGraphs.AbstractGraph, P_inits::Vector{Any}, C::AbstractFloat, annealing_schedule::Function, annealing_schedule_name::String, steps_per_temp::Integer, k_max::Integer, N_runs::Integer)
    energies = []
    P_finals = []
    for i in 1:N_runs
        P, en = sim_anneal(g, P_inits[i], C, annealing_schedule, steps_per_temp, k_max)
        push!(energies, en)
        push!(P_finals, P)
    end
    N_vertices = length(P_inits[1])
    JLD.save(filepath, "energies",energies, "P_inits",P_inits, "P_finals",P_finals, "N_vertices",N_vertices, "Grid",g, "annealing_schedule",annealing_schedule_name, "steps_per_temp",steps_per_temp, "C",C , "k_max",k_max, "N_runs",N_runs)
end

using Distributed
""" For parallel computing on high performance clusters.
"""
function parallel_multiple_sim_anneal(filepath::String, g::LightGraphs.AbstractGraph, P_inits::Vector{Any}, C::AbstractFloat, annealing_schedule::Function, annealing_schedule_name::String, steps_per_temp::Integer, k_max::Integer, N_runs::Integer)
    Data = @distributed vcat for i in 1:N_runs
        sim_anneal(g, P_inits[i], C, annealing_schedule, steps_per_temp, k_max)
        #P, en = sim_anneal(g, P_inits[i], C, annealing_schedule, steps_per_temp, k_max)
        #push!(energies, en)
        #push!(P_finals, P)
    end
    N_vertices = length(P_inits[1])
    P_finals = first.(Data)
    energies = last.(Data)
    JLD.save(filepath, "energies",energies, "P_inits",P_inits, "P_finals",P_finals, "N_vertices",N_vertices, "Grid",g, "annealing_schedule",annealing_schedule_name, "steps_per_temp",steps_per_temp, "C",C , "k_max",k_max, "N_runs",N_runs)
end


### several functions for SA
function probability(ΔE::AbstractFloat, T::AbstractFloat) # probability function depends on ΔE and on temperature function
    exp( - ΔE / T) # Boltzmann's constant k_B is set to one
end

################################################################################
################### temperature schedules for annealing  #######################
################################################################################

# arbitrary temperature function that decreases to zero and calculates temperature dependent of step
function temp1(k::Integer, steps_per_temp::Integer)
    1. / (1 + 0.0001 * floor(k/steps_per_temp))
end

function temp2(k::Integer, steps_per_temp::Integer)
    1. / (1 + 0.0007 * floor(k/steps_per_temp))
end

function temp3(k::Integer, steps_per_temp::Integer)
    1. / (1 + 0.0015 * floor(k/steps_per_temp))
end

function temp4(k::Integer, steps_per_temp::Integer)
    1. / (1 + 0.0025 * floor(k/steps_per_temp))
end

function temp_ex1_a(k::Integer, steps_per_temp::Integer)
    0.99 ^ (floor(k/steps_per_temp)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex1_b(k::Integer, steps_per_temp::Integer)
    0.99 ^ (floor(k/steps_per_temp)) + 0.25 # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex1_c(k::Integer, steps_per_temp::Integer)
    0.99 ^ (floor(k/steps_per_temp)) + 0.5 # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex2(k::Integer, steps_per_temp::Integer)
    0.999 ^ (floor(k/steps_per_temp)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end


################################################################################
######################### DEPRECATED FUNCTIONS##################################
################################################################################

# # core function for simulated annealing
# function sim_anneal(g::LightGraphs.AbstractGraph, P_init::Array{Float64,1}, C::AbstractFloat, annealing_schedule::Function, steps_per_temp::Integer, k_max::Integer) # k_max: setting number of computation steps
#     # given an initial configuration P sim_anneal() tendentially finds a more stable configuration
#
#     P = copy(P_init)
#     en = [ ]
#     for k in 0:k_max - 1
#         T = annealing_schedule(k,steps_per_temp) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
#         P_old = copy(P) # for calculating the energy of "old" configuration
#         P = stable_swapped_config!(g, P, C)
#         energy_old = energy(g, P_old, C)[2] # by [2] only the second value of tuple is returned (G_av)
#         energy_new = energy(g, P, C)[2] # by [2] only the second value of tuple is returned (G_av)
#         ΔE = energy_new - energy_old
#         #### performance: let energy() calculate G_av only? Nope, G is only saving
#         # to an array and is helpful for understanding the algorithm.
#
#         if ΔE <= 0 # man könnte conditional auch umdrehen: if (ΔE <= 0 AND probability(ΔE, T) < rand())
#                                                                  # P = P_old
#             P
#         elseif probability(ΔE, T) > rand() # rand() gives random number element of [0,1]
#             P
#         else
#             P = P_old
#         end
#         push!(en, energy_old)
#     end
#     P, en
# end
