struct Locality
    loc_1step_init
    loc_1step_final
    loc_1sept_0_init
    loc_1sept_0_final
end

struct Nr_gen_con
    gen_gen
    con_con
    gen_con
end

""" Function that takes .jld as input containing the value of the energy for each
    iteration step of simulated annealing and executes multiple postprocessing
    functions on that energy values. The resulting postprocessing data is save as
    another .jld
"""
function postprocess_sim_anneal(filepath_in, filepath_out, T)
    Data_loaded = JLD.load(filepath_in)
    energies = Data_loaded["energies"]
    P_init = Data_loaded["P_init"]
    P_finals = Data_loaded["P_finals"]
    N_vertices = Data_loaded["N_vertices"]
    g = Data_loaded["Grid"]
    ann_sched = Data_loaded["annealing_schedule"]
    steps_per_temp = Data_loaded["steps_per_temp"]
    C = Data_loaded["C"]
    k_max = Data_loaded["k_max"]
    N_runs = Data_loaded["N_runs"]

    energy_init = []
    energy_final = []
    N_T_init = []
    N_T_final = []
    locality = []
    #locality_0 = []
    nr_gen_con_init = []
    nr_gen_con_final = []

    for i in 1:N_runs

        # calculating observables
        # energy of initial and final configutations
        push!(energy_init, energies[i][1])
        push!(energy_final, energies[i][k_max])

        # number of flows above a certain threshold for initial and final configutations
        push!(N_T_init, flows_above_thres(T, P_init, g))
        push!(N_T_init, flows_above_thres(T, P_finals[i], g))

        # locality
        # locality_init_final = loc_1step(g, P_init, C), loc_1step(g, P_finals[i], C)
        # push!(locality, locality_init_final)
        # locality_0_init_final = loc_1step_0(g, P_init, C), loc_1step_0(g, P_finals[i], C)
        # push!(locality_0, locality_0_init_final)
        locality_single_run = Locality(loc_1step(g, P_init, C),loc_1step(g, P_finals[i], C),loc_1step_0(g, P_init, C),loc_1step_0(g, P_finals[i], C))
        push!(locality, locality_single_run)

        # nr_gen_con
        gen_gen_init, con_con_init, gen_con_init = nr_gen_con(g,P_init)
        nr_gen_con_single_run_init = Nr_gen_con(gen_gen_init, con_con_init, gen_con_init)
        push!(nr_gen_con_init,nr_gen_con_single_run_init)
        gen_gen_final, con_con_final, gen_con_final = nr_gen_con(g,P_finals[i])
        nr_gen_con_single_run_final = Nr_gen_con(gen_gen_final, con_con_final, gen_con_final)
        push!(nr_gen_con_final,nr_gen_con_single_run_final)

    end
    JLD.save(filepath_out, "energies",energies, "P_init",P_init, "P_finals",P_finals, "N_vertices",N_vertices, "Grid",g,
        "annealing_schedule",ann_sched, "steps_per_temp",steps_per_temp, "C",C , "k_max",k_max, "N_runs",N_runs,
        "energy_init",energy_init, "energy_final",energy_final, "N_T_init",N_T_init, "N_T_final",N_T_final,
        "locality",locality, "nr_gen_con_init",nr_gen_con_init, "nr_gen_con_final",nr_gen_con_final)

        #, "locality_0",locality_0)
end


# histograms
function flow_single_sample(SData, i, j) # i: sample number, j = 1: random, j = 4: optimized
    N_side = SData["N_side"]
    g = gen_square_grid(N_side)
    Data = SData["Data"]
    P = Data[i][j]
    F = flow(g, P)
    x = abs.(F)
end

function flows_N_runs(SData, j) # i: sample number, j = 1: random grids, j = 4: optimised grids
    N_side = SData["N_side"]
    g = gen_square_grid(N_side)
    Data = SData["Data"]
    All_flows = [ ]
    N_runs = length(Data)
    P = Data[1][j]
    F = flow(g, P)
    All_flows = copy(F)
    for i in 2:N_runs
        P = Data[i][j]
        F = flow(g, P)
        append!(All_flows, F)
    end
    All_flows
    x = abs.(All_flows)
end

function high_gc_low_Gav(SData, gen_con, G_av_final)
    Data = SData["Data"]
    N_runs = length(Data)
    configs_high_gc_low_Gav = [ ]
    for i in 1:N_runs
        if Data[i][6][3] > gen_con && Data[i][5][2] < G_av_final # 23 is mean G_av plus STD
            configs_high_gc_low_Gav = push!(configs_high_gc_low_Gav, Data[i])
        end
    end
    configs_high_gc_low_Gav
end

function flows_high_gc_low_Gav(Data, N_side, j) # i: sample number, j = 1: random grids, j = 4: optimised grids
    g = gen_square_grid(N_side)
    All_flows = [ ]
    N_runs = length(Data)
    P = Data[1][j]
    F = flow(g, P)
    All_flows = copy(F)
    for i in 2:N_runs
        P = Data[i][j]
        F = flow(g, P)
        append!(All_flows, F)
    end
    All_flows
    x = abs.(All_flows)
end

function weights_mean_err(SData, j, nbins) # j = 1: random grids, j = 4: optimised grids
    weights = [ ]
    N_bars = nbins - 1
    for i in 1:N_bars
        push!(weights, [ ])
    end
    Data = SData["Data"]
    N_runs = length(Data)
    bins = range(0.0, 1.0, length = nbins)
    for i in 1:N_runs
        flows = flow_single_sample(SData, i, j)
        h = fit(Histogram, flows, bins)
        h = normalize(h, mode=:probability)
        weightvalues = h.weights
        for i in 1:N_bars
            append!(weights[i], weightvalues[i])
        end
    end
    weights_av = [ ]
    weights_err = [ ]
    for i in 1:N_bars
        append!(weights_av, mean(weights[i]))
        append!(weights_err, 1 / sqrt(N_runs) * Statistics.std(weights[i]))
    end
    weights_av, weights_err
end

""" Calculates number of edges gen-gen, con-con, gen-con.
"""
function nr_gen_con(g, P)
    N_edges = ne(g)
    B = Array(incidence_matrix(g, oriented=true))
    gen_gen = 0 # generates empty vector of type 'Any'
    con_con = 0
    gen_con = 0

    # go through all edges
    for i in 1:N_edges
        # for each edge find the two indizes of the two vertices that are connected by that edge
        index_P_first = findfirst(isodd, B[:, i])
        index_P_second = findlast(isodd, B[:, i])
        # read out value for each vertex and check if this edge is connected by gen_gen, con_con or gen_con
        if P[index_P_first] == 1 && P[index_P_second] == 1
            gen_gen = gen_gen + 1
        elseif P[index_P_first] == -1 && P[index_P_second] == -1
            con_con = con_con + 1
        else gen_con = gen_con + 1
        end
    end
    gen_gen, con_con, gen_con
end


function nr_gen_con_av(SData, col)
    Data = SData["Data"]
    N_runs = SData["N_runs"]
    gen_gen_av = mean(collect_data2(Data, col, 1))
    con_con_av = mean(collect_data2(Data, col, 2))
    gen_con_av = mean(collect_data2(Data, col, 3))
    gen_gen_std = 1 / sqrt(N_runs) * std(collect_data2(Data, col, 1))
    con_con_std = 1 / sqrt(N_runs) * std(collect_data2(Data, col, 2))
    gen_con_std = 1 / sqrt(N_runs) * std(collect_data2(Data, col, 3))
    gen_gen_av, gen_gen_std, con_con_av, con_con_std, gen_con_av, gen_con_std
end

function nr_gen_con_av_diff(SData)
    Data = SData["Data"]
    N_runs = SData["N_runs"]
    gen_gen_av = mean(collect_data2(Data, 3, 1) - collect_data2(Data, 6, 1))
    con_con_av = mean(collect_data2(Data, 3, 2) - collect_data2(Data, 6, 2))
    gen_con_av = mean(collect_data2(Data, 3, 3) - collect_data2(Data, 6, 3))
    gen_gen_std = sqrt((1 / sqrt(N_runs) * std(collect_data2(Data, 3, 1)))^2 + (1 / sqrt(N_runs) * std(collect_data2(Data, 6, 1)))^2) # the measurement uncertainties sum up
    con_con_std = sqrt((1 / sqrt(N_runs) * std(collect_data2(Data, 3, 2)))^2 + (1 / sqrt(N_runs) * std(collect_data2(Data, 6, 2)))^2)
    gen_con_std = sqrt((1 / sqrt(N_runs) * std(collect_data2(Data, 3, 3)))^2 + (1 / sqrt(N_runs) * std(collect_data2(Data, 6, 3)))^2)
    gen_gen_av, gen_gen_std, con_con_av, con_con_std, gen_con_av, gen_con_std
end


# locality functions
""" Removes all edges whose flow is bigger than C -> calculates flow again and so on
    until all flows are smaller than C, saves intermediate graphs.
"""
function cascading_steps!(g, P, C)
    F = flow(g, P)
    cascading_steps = [ ]
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
        g_new = copy(g)
        push!(cascading_steps, g_new)
        F = flow(g, P)
    end
    cascading_steps
end

""" Removes all edges whose flow is bigger than C, outputs vertices of those edges.
"""
function cascade_1step!(g, P, C) #
    F = flow(g, P)
    if maximum(abs.(F)) > C
        B = Array(incidence_matrix(g, oriented=true))
        n = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
        # it is ne(g) = size(B)[2], where ne() give number of edges
        vertices_of_failed_edges = [ ]
        for i in 1:n # in this loop all edges that carry flow bigger than certain value/ threshold are being removed
            if abs(F[i]) > C
                rem_edge!(g, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
                # B[:, i] gives elements of i-th column of B, isodd accesses only the odd elements of these element
                # findfirst()/ findlast() access first/ last odd element.
                # incidence matrix B comprises only zeros (even elements) and one 1 and one -1 for the two vertices that are connected
                # to respective edge
                push!(vertices_of_failed_edges, findfirst(isodd, B[:, i]), findlast(isodd, B[:, i]))
            end
        end
    else vertices_of_failed_edges = [ ]
    end
    vertices_of_failed_edges
end

""" Calculates distance from initial line failure to closest secondary failure. This
    distance is calculated for every edge and then the average of the distances of
    all edges is calculated.
"""
function loc_1step(g_init, P, C)
    g = copy(g_init)
    B = Array(incidence_matrix(g, oriented=true))
    n = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns  i.e. the number of edges
    linefailure_indizes = collect(1:n) # collect() collects all values in the range 1:n in an array, here all edges are numberd in an array
    #linefailure_indizes = shuffle!(linefailure_indizes) # positions of edges in linefailure_indizes is randomly permuted by shuffle!()

    min_failure_distances = [ ]
    for i in 1:n
        g = copy(g_init)
        B = Array(incidence_matrix(g, oriented=true))
        first = findfirst(isodd, B[:, i]) # finds first element that is 1
        last = findlast(isodd, B[:, i]) # finds last element that is 1
        rem_edge!(g, first, last)
        vertices_initial_failed_edge = [first, last]
        vertices_of_failed_edges = cascade_1step!(g, P, C)
        if isempty(vertices_of_failed_edges)
            continue
        end
        failure_distances = [ ]
        for i in 1:length(vertices_initial_failed_edge)
            for j in 1:length(vertices_of_failed_edges)
                push!(failure_distances, length(a_star(g,vertices_initial_failed_edge[i],vertices_of_failed_edges[j])))
            end
        end
        push!(min_failure_distances, minimum(failure_distances))
    end
    if isempty(min_failure_distances)
        min_failure_distances = [0]
    end
    min_failure_distance_av = mean(min_failure_distances)
    min_failure_distances, min_failure_distance_av
end

function loc_1step_0(g_init, P, C) # only differing to loc_1step!() by counting the distance
# as zero in case no other edge failure is caused by initial edge removal
    g = copy(g_init)
    B = Array(incidence_matrix(g, oriented=true))
    n = size(B)[2] # size() gives array containing number of rows and columns of incidence matrix b, [2] accesses number of columns
    linefailure_indizes = collect(1:n) # collect() collects all values in the range 1:m in an array, here all edges are numberd in an array
    #linefailure_indizes = shuffle!(linefailure_indizes) # positions of edges in linefailure_indizes is randomly permuted by shuffle!()

    min_failure_distances = [ ]
    for i in 1:n
        g = copy(g_init)
        B = Array(incidence_matrix(g, oriented=true))
        first = findfirst(isodd, B[:, i])
        last = findlast(isodd, B[:, i])
        rem_edge!(g, first, last)
        vertices_initial_failed_edge = [first, last]
        vertices_of_failed_edges = cascade_1step!(g, P, C)
        failure_distances = [ ]
        if isempty(vertices_of_failed_edges)
            push!(failure_distances, 0)
        end

        for i in 1:length(vertices_initial_failed_edge)
            for j in 1:length(vertices_of_failed_edges)
                push!(failure_distances, length(a_star(g,vertices_initial_failed_edge[i],vertices_of_failed_edges[j])))
            end
        end
        push!(min_failure_distances, minimum(failure_distances))
    end
    min_failure_distance_av = mean(min_failure_distances)
    min_failure_distances, min_failure_distance_av
end

function locality(SData)
    Data = SData["Data"]
    N_side = SData["N_side"]
    C = SData["C"]
    N_runs = SData["N_runs"]
    configs_init= collect_data(Data, 1)
    configs_final = collect_data(Data, 4)
    Data_init = [ ]
    Data_final = [ ]
    Data_init0 = [ ]
    Data_final0 = [ ]
    N_runs = length(Data)
    for i in 1:N_runs
        x = loc_1step!(configs_init[i], C, N_side)
        push!(Data_init, x)
        y = loc_1step!(configs_final[i], C, N_side)
        push!(Data_final, y)
        v = loc_1step_0!(configs_init[i], C, N_side)
        push!(Data_init0, v)
        w = loc_1step_0!(configs_final[i], C, N_side)
        push!(Data_final0, w)
    end
    locality_av = mean(collect_data(Data_final, 2) - collect_data(Data_init, 2))
    locality_std = sqrt((1 / sqrt(N_runs) * std(collect_data(Data_final, 2)))^2 + (1 / sqrt(N_runs) * std(collect_data(Data_init, 2)))^2)
    locality_av0 = mean(collect_data(Data_final0, 2) - collect_data(Data_init0, 2))
    locality_std0 = sqrt((1 / sqrt(N_runs) * std(collect_data(Data_final0, 2)))^2 + (1 / sqrt(N_runs) * std(collect_data(Data_init0, 2)))^2)
    #locality_std0 = 1 / sqrt(N_runs) * (std(collect_data(Data_final0, 2)) + std(collect_data(Data_init0, 2)))
    A = locality_av, locality_std, Data_init, Data_final, mean(collect_data(Data_init, 2)), 1 / sqrt(N_runs) * std(collect_data(Data_init, 2)), mean(collect_data(Data_final, 2)), 1 / sqrt(N_runs) * std(collect_data(Data_final, 2))
    B = locality_av0, locality_std0, Data_init0, Data_final0, mean(collect_data(Data_init0, 2)), 1 / sqrt(N_runs) * std(collect_data(Data_final0, 2)), mean(collect_data(Data_final0, 2)), 1 / sqrt(N_runs) * std(collect_data(Data_init0, 2))
    A, B
end


################################################################################
######################### DEPRECATED FUNCTIONS##################################
################################################################################

# # measure  number of edges gen-gen, con-con, gen-con 04.05.2020
# # This function holds only for square grids!
# function nr_gen_con_old(P, N_side)
#     N_vertices = N_side * N_side
#     gen_gen = 0 # generates empty vector of type 'Any'
#     con_con = 0
#     gen_con = 0
#     for i in 1:N_side:N_vertices - N_side + 1
#         for i in i:i + N_side - 2
#             if P[i] == 1 && P[i + 1] == 1
#                 gen_gen = gen_gen + 1
#             elseif P[i] == -1 && P[i + 1] == -1
#                 con_con = con_con + 1
#             else gen_con = gen_con + 1
#             end
#         end
#     end
#     for i in 1:N_side
#         for i in i:N_side:N_vertices - (2 * N_side) + i
#             if P[i] == 1 && P[i + N_side] == 1
#                 gen_gen = gen_gen + 1
#             elseif P[i] == -1 && P[i + N_side] == -1
#                 con_con = con_con + 1
#             else gen_con = gen_con + 1
#             end
#         end
#     end
#     gen_gen, con_con, gen_con
# end
