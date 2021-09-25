#= This file contains functions that calculate observerbles out of the "raw" simulation
data. =#

struct Locality
    loc_1step_init
    loc_1step_final
    loc_1step_0_init
    loc_1step_0_final
end

struct Nr_gen_con
    gen_gen
    con_con
    gen_con
end

""" Function that takes .jld-file as input containing the value of the energy for
    each iteration step of simulated annealing and executes multiple postprocessing
    functions. The resulting postprocessing data is saved as another .jl-file.
"""
function postprocess_sim_anneal(filepath_in, filepath_out, T)
    Data_loaded = JLD.load(filepath_in)
    energies = Data_loaded["energies"]
    P_inits = Data_loaded["P_inits"]
    P_finals = Data_loaded["P_finals"]
    N_vertices = Data_loaded["N_vertices"]
    g = Data_loaded["Grid"]
    ann_sched = Data_loaded["annealing_schedule"]
    steps_per_temp = Data_loaded["steps_per_temp"]
    C = Data_loaded["C"]
    k_max = Data_loaded["k_max"]
    N_runs = Data_loaded["N_runs"]

    energy_init = []; energy_final = []
    N_T_init = []; N_T_final = []
    loc_1step_init = []; loc_1step_final = []
    loc_1step_0_init = []; loc_1step_0_final = []
    gen_gen_init = []; con_con_init = []; gen_con_init = []
    gen_gen_final = []; con_con_final = []; gen_con_final = []

    for i in 1:N_runs

        # calculating observables
        # energy of initial and final configutations
        push!(energy_init, energies[i][1])
        push!(energy_final, energies[i][k_max])

        # number of flows above a certain threshold for initial and final configutations
        push!(N_T_init, flows_above_thres(T, P_inits[i], g))
        push!(N_T_final, flows_above_thres(T, P_finals[i], g))

        # locality
        loc_1step_init_single_run = loc_1step(g, P_inits[i], C)
        loc_1step_final_single_run = loc_1step(g, P_finals[i], C)
        loc_1step_0_init_single_run = loc_1step_0(g, P_inits[i], C)
        loc_1step_0_final_single_run = loc_1step_0(g, P_finals[i], C)
        push!(loc_1step_init, loc_1step_init_single_run)
        push!(loc_1step_final, loc_1step_final_single_run)
        push!(loc_1step_0_init, loc_1step_0_init_single_run)
        push!(loc_1step_0_final, loc_1step_0_final_single_run)

        # nr_gen_con
        gen_gen_single_run_init, con_con_single_run_init, gen_con_single_run_init = nr_gen_con(g,P_inits[i])
        push!(gen_gen_init,gen_gen_single_run_init); push!(con_con_init,con_con_single_run_init); push!(gen_con_init,gen_con_single_run_init);
        gen_gen_single_run_final, con_con_single_run_final, gen_con_single_run_final = nr_gen_con(g,P_finals[i])
        push!(gen_gen_final,gen_gen_single_run_final); push!(con_con_final,con_con_single_run_final); push!(gen_con_final,gen_con_single_run_final);
    end

    nr_gen_con_init = Nr_gen_con(gen_gen_init, con_con_init, gen_con_init)
    nr_gen_con_final = Nr_gen_con(gen_gen_final, con_con_final, gen_con_final)
    locality = Locality(loc_1step_init,loc_1step_final,loc_1step_0_init,loc_1step_0_final)

    JLD.save(filepath_out, "energies",energies, "P_inits",P_inits, "P_finals",P_finals, "N_vertices",N_vertices, "Grid",g,
        "annealing_schedule",ann_sched, "steps_per_temp",steps_per_temp, "C",C , "k_max",k_max, "N_runs",N_runs,
        "energy_init",energy_init, "energy_final",energy_final, "N_T_init",N_T_init, "N_T_final",N_T_final,
        "locality",locality, "nr_gen_con_init",nr_gen_con_init, "nr_gen_con_final",nr_gen_con_final)
end


### results SA: calculating Gav_av and STD_Gav (averaged over all runs)
function Gav_av_STD_Gav(Data_loaded)
    N_runs = Data_loaded["N_runs"]
    # random grids
    Gav_av_init = round.(mean(Data_loaded["energy_init"]); digits = 2)
    STD_Gav_init = round.(1 / sqrt(N_runs) * std(Data_loaded["energy_init"]; corrected=true); digits = 2) # corrected=true is default of std(), so it could be omitted
    # minimized G_av
    Gav_av_final = round.(mean(Data_loaded["energy_final"]); digits = 2)
    STD_Gav_final = round.(1 / sqrt(N_runs) * std(Data_loaded["energy_final"]; corrected=true); digits = 2)
    Gav_av_init, STD_Gav_init, Gav_av_final, STD_Gav_final
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


function nr_gen_con_av(Data_loaded, P_rand_opt)
    N_runs = Data_loaded["N_runs"]
    gen_gen_av = mean(Data_loaded[P_rand_opt].gen_gen)
    con_con_av = mean(Data_loaded[P_rand_opt].con_con)
    gen_con_av = mean(Data_loaded[P_rand_opt].gen_con)
    gen_gen_std = 1 / sqrt(N_runs) * std(Data_loaded[P_rand_opt].gen_gen)
    con_con_std = 1 / sqrt(N_runs) * std(Data_loaded[P_rand_opt].con_con)
    gen_con_std = 1 / sqrt(N_runs) * std(Data_loaded[P_rand_opt].gen_con)
    gen_gen_av, gen_gen_std, con_con_av, con_con_std, gen_con_av, gen_con_std
end

function nr_gen_con_av_diff(Data_loaded)
    N_runs = Data_loaded["N_runs"]
    gen_gen_av = mean(Data_loaded["nr_gen_con_init"].gen_gen - Data_loaded["nr_gen_con_final"].gen_gen)
    con_con_av = mean(Data_loaded["nr_gen_con_init"].con_con - Data_loaded["nr_gen_con_final"].con_con)
    gen_con_av = mean(Data_loaded["nr_gen_con_init"].gen_con - Data_loaded["nr_gen_con_final"].gen_con)
    gen_gen_std = sqrt((1 / sqrt(N_runs) * std(Data_loaded["nr_gen_con_init"].gen_gen))^2 + (1 / sqrt(N_runs) * std(Data_loaded["nr_gen_con_final"].gen_gen))^2) # the measurement uncertainties sum up
    con_con_std = sqrt((1 / sqrt(N_runs) * std(Data_loaded["nr_gen_con_init"].con_con))^2 + (1 / sqrt(N_runs) * std(Data_loaded["nr_gen_con_final"].con_con))^2)
    gen_con_std = sqrt((1 / sqrt(N_runs) * std(Data_loaded["nr_gen_con_init"].gen_con))^2 + (1 / sqrt(N_runs) * std(Data_loaded["nr_gen_con_final"].gen_con))^2)
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

function locality(Data_loaded)
    C = Data_loaded["C"]
    N_runs = Data_loaded["N_runs"]
    Data_init = [ ]
    Data_final = [ ]
    Data_init0 = [ ]
    Data_final0 = [ ]
    for i in 1:N_runs
        x = Data_loaded["locality"].loc_1step_init[i][2]
        push!(Data_init, x)
        y = Data_loaded["locality"].loc_1step_final[i][2]
        push!(Data_final, y)
        v = Data_loaded["locality"].loc_1step_0_init[i][2]
        push!(Data_init0, v)
        w = Data_loaded["locality"].loc_1step_0_final[i][2]
        push!(Data_final0, w)
    end
    # calculating difference between random and optimized configurations
    locality_av_diff = mean(Data_final - Data_init)
    locality_std_diff = sqrt((1 / sqrt(N_runs) * std(Data_final))^2 + (1 / sqrt(N_runs) * std(Data_init))^2)
    locality_av0_diff = mean(Data_final0 - Data_init0)
    locality_std0_diff = sqrt((1 / sqrt(N_runs) * std(Data_final0))^2 + (1 / sqrt(N_runs) * std(Data_init0))^2)
    #locality_std0 = 1 / sqrt(N_runs) * (std(collect_data(Data_final0, 2)) + std(collect_data(Data_init0, 2)))
    A = locality_av_diff, locality_std_diff, Data_init, Data_final, mean(Data_init), 1 / sqrt(N_runs) * std(Data_init), mean(Data_final), 1 / sqrt(N_runs) * std(Data_final)
    B = locality_av0_diff, locality_std0_diff, Data_init0, Data_final0, mean(Data_init0), 1 / sqrt(N_runs) * std(Data_final0), mean(Data_final0), 1 / sqrt(N_runs) * std(Data_init0)
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
