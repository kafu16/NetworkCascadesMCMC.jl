################################################################################
##### grid operations (funtions that change or refer to  grid topology) ########
################################################################################


### failure/removal of edge that causes cascade
function linefailure!(g, i) # edge i is removed
    B = Array(incidence_matrix(g, oriented=true))
    rem_edge!(g, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
    g
end


### cascade
function cascade!(g, P, C) # removes all edges whose flow is bigger than C -> calculates flow again and so on until all flows are smaller than C
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
function biggest_component(g) # gives biggest connected component of graph g
    maximum(length.(LightGraphs.connected_components(g))) # connected_components() returns vector whose components contain connected vertices
    # so the length of a component is the size of a connected graph
    # so here the biggest connected graph out of all vertices is chosen
end


################################################################################
############################# energy function(s) ###############################
################################################################################
using Statistics

function energy(g_init, P, C, N_side) # calculates energy of step k, C: threshold that marks line failure,
    g = copy(g_init)
    B = Array(incidence_matrix(g, oriented=true))
    m = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
    linefailure_indizes = collect(1:m) # collect() collects all values in the range 1:m in an array, here all edges are numberd in an array
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

    G = zeros(m)
    for i in 1:m # for loop for randomly chosen N_vertices linefailures
        g = linefailure!(g, linefailure_indizes[i])
        g = cascade!(g, P, C)
        #global G[i] = ne(g)
        global G[i] = biggest_component(g) # G: size of biggest connected component, global lets the values of G saved after each run of loop,
        # otherwise G would be overwritten in each run of loop
        g = copy(g_init)
    end

    G_av = mean(G)
    G, G_av # this way two values in a tuple are returned by a function
end


################################################################################
########################## Monte Carlo step functions ##########################
################################################################################

### swap of configurations
function swap!(P, N_side) # swaps randomly chosen generator-consumer pair
    N_vertices = length(P)
    A = rand(1:N_vertices,1)
    B = rand(1:N_vertices,1)
    # println("A initial", A) # following println-lines only to check what while-loop does, not necessary for swap()
    # println("B initial", B)
    while A == B || P[A] == P[B]
    # A == B avoids that entity is exchanged by itself, P[A] = P[B] avoids consumer-consumer and generator-generator-swap,
    # while-loop is run if condition on either left or right hand side of || is fulfilled
        # println("A while condition", A, "P[A] while condition", P[A])
        # println("B while condition", B, "P[B] while condition", P[B])
        A = rand(1:N_vertices,1) # if one writes "!" after function, behaviour of function changes see https://docs.julialang.org/en/v1/stdlib/Random/
        B = rand(1:N_vertices,1) #### ToDo warum Fehlermeldung, wenn man A und B als global assigned?
        # println("A in while", A)
        # println("B in while", B)
    end
    # println("A after while", A, "P[A] after while", P[A])
    # println("B after while", B, "P[B] after while", P[B])
    P[A], P[B] = P[B], P[A] # swaps neighbors
    # println(P)
    P
end


# generation of stable swapped configuration
function stable_swapped_config!(g, P, C, N_side)
# to avoid calculation steps, the input configuration should be stable as the stable_swapped_config() only permutes
# one generator-consumer pair. so having an unstable configuration as input will probably take more steps than
# first generating a stable configuration by gen_stable_config() and then apply stable_swapped_config()
#### ToDo build: return error when stable config is not possible due to a too low C
    P_stable_old = copy(P)
    # by using P_stable_old it is made sure that a the swapped configuration differs only by one permutation
    # otherwise by the following while-loop in each iteration of the loop the newly generated configuration would be swapped
    # again to obtain a stable configuration. This can not be done by P_stable_old = P because of variable assignment. In the latter case
    # P_stable_old would equal P all the time, even when the value of P is changed
    P = swap!(P, N_side)
    F = flow(g, P)
    #### ToDo: include max_iterations as variable parameter?
    max_iterations = 10000 # hardcoded number of iteration after that loop is excited
    global i = 1
    while maximum(abs.(F)) > C
        P = swap!(P_stable_old, N_side)
        F = flow(g, P)
        #### ToDo
        if i >= max_iterations
            error("ERROR: Maximum number of iterations for finding stable configuration reached.")
        end
        global i += 1
    end
    P
end


### evaluation of energy one random configuration
N_removals = 0
# """
# function for evaluating the biggest component without evaluating several other functions
# """
function eval_one_random_config(N_side, C)
    g = gen_square_grid(N_side)
    P = gen_stable_config(g, N_side, C)
    energy(g, P, C, N_side)
end


################################################################################
############### implementation of Simulated Annealing (SA) #####################
################################################################################

### key code of SA

#### ToDo: Option N_removals entfernen

# core function for simulated annealing
function sim_anneal!(g, P, C, N_side, k_max) # k_max: setting number of computation steps
    # given an initial configuration P sim_anneal!() tendentially finds a more stable configuration

    en = [ ]
    for k in 0:k_max - 1
        T = temperature(k) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
        P_old = copy(P) # for calculating the energy of "old" configuration
        P = stable_swapped_config!(g, P, C, N_side)
        energy_old = energy(g, P_old, C, N_side)[2] # by [2] only the second value of tuple is returned (G_av)
        energy_new = energy(g, P, C, N_side)[2] # by [2] only the second value of tuple is returned (G_av)
        ΔE = energy_new - energy_old
        #### performance: let energy() calculate G_av only? Nope, G is only saving
        # to an array and is helpful for understanding the algorithm.

        if ΔE <= 0 # man könnte conditional auch umdrehen: if (ΔE <= 0 AND probability(ΔE, T) < rand())
                                                                 # P = P_old
            P
        elseif probability(ΔE, T) > rand() # rand() gives random number element of [0,1]
            P
        else
            P = P_old
        end
        push!(en, energy_old)
    end
    P, en
end

# applies SA on random square grid
#### ToDo evtl. diesen code direkt in sim_anneal() einbauen. Bin mir unsicher, ob das sinnvoll ist...
function eval_sim_anneal!(N_side, C, T, N_removals = 0, k_max = 10)
    g = gen_square_grid(N_side)
    P = gen_stable_config(g, N_side, C) # to avoid iteration steps it is important to start with a stable configurations see comment at stable_swapped_config!()
    P_initial = copy(P)
    P, en = sim_anneal!(g, P, C, N_side, 0, k_max)
    g = gen_square_grid(N_side)
    energy_initial = energy(g, P_initial, C, N_side)
    g = gen_square_grid(N_side)
    N_T = flows_above_thres(T, P, g)
    energy_final = energy(g, P, C, N_side)
    P_initial, energy_initial, P, energy_final, N_T, en
end

### several functions for SA

function probability(ΔE, T) # probability function depends on ΔE and on temperature function
    exp( - ΔE / T) # Boltzmann's constant k_B is set to one
end

################################################################################
################### temperature schedules for annealing  #######################
################################################################################

# arbitrary temperature function that decreases to zero and calculates temperature dependent of step
function temperature(k)
    0.99 ^ k
    #0.999 ^ (floor(k/4))
    #1. / (1 + 0.0001 * floor(k/4))
end

# arbitrary temperature function that decreases to zero and calculates temperature dependent of step
function temp1(k)
    1. / (1 + 0.0001 * floor(k/4))
end

function temp2(k)
    1. / (1 + 0.0007 * floor(k/4))
end

function temp3(k)
    1. / (1 + 0.0015 * floor(k/4))
end

function temp4(k)
    1. / (1 + 0.0025 * floor(k/4))
end

function temp5(k)
    10. / (1 + 0.02 * k)
end

function temp_ex1(k)
    0.999 ^ (floor(k/4)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex2(k)
    0.9995 ^ (floor(k/4)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex3(k)
    0.99975 ^ (floor(k/4)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex4(k)
    0.9999 ^ (floor(k/4)) # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end

function temp_ex5(k)
    0.99 ^ k # floor(x) returns the nearest integral value of the same type as x that is less than or equal to x
end
