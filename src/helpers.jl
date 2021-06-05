#= this file contains helper functions that don't fit in any other file of
NetworkCascadesMCMC=#

# data collection
function collect_data_SA_runs_var_ann_shed(N_runs, N_side, C, T, annealing_schedule, k_max)
    Data = []
    for i in 1:N_runs
        g = gen_square_grid(N_side)
        P = gen_stable_config(g, N_side, C) # to avoid iteration steps it is important to start with a stable configurations see comment at stable_swapped_config!()
        P_initial = copy(P)
        en = [ ]
        N_removals = 0
        for k in 0:k_max - 1
            Temp = annealing_schedule(k)
            P_old = copy(P) # for calculating the energy of "old" configuration
            P = stable_swapped_config!(g, P, C)
            energy_old = energy(g, P_old, C)[2] # by [2] only the second value of tuple is returned (G_av)
            energy_new = energy(g, P, C)[2] # by [2] only the second value of tuple is returned (G_av)
            ΔE = energy_new - energy_old
            #### performance: let energy() calculate G_av only

            if ΔE <= 0 # man könnte conditional auch umdrehen: if (ΔE <= 0 AND probability(ΔE, T) < rand())
                                                                     # P = P_old
                P
            elseif probability(ΔE, Temp) > rand() # rand() gives random number element of [0,1]
                P
            else
                P = P_old
            end
            g = gen_square_grid(N_side)
            push!(en, energy_old)
        end
        energy_initial = energy(g, P_initial, C)
        N_T = flows_above_thres(T, P_initial, g), flows_above_thres(T, P, g)
        locality_init = loc_1step!(P_initial, C, N_side), loc_1step_0!(P_initial, C, N_side)
        locality_final = loc_1step!(P, C, N_side), loc_1step_0!(P, C, N_side)
        energy_final = energy(g, P, C)
        SA_extremal = P_initial, energy_initial, nr_gen_con(P_initial, N_side), P, energy_final, nr_gen_con(P, N_side), N_T, en, locality_init, locality_final
        push!(Data, SA_extremal)
    end
    Data
end


# function collect_data_SA_runs_var_ann_shed(N_runs, N_side, C, T, annealing_schedule, k_max)
#     Data = []
#     for i in 1:N_runs
#         g = gen_square_grid(N_side)
#         P = gen_stable_config(g, N_side, C) # to avoid iteration steps it is important to start with a stable configurations see comment at stable_swapped_config!()
#         P_initial = copy(P)
#         P, en = sim_anneal!(g, P, C, 0, k_max)
#         energy_initial = energy(g, P_initial, C)
#         N_T = flows_above_thres(T, P_initial, g), flows_above_thres(T, P, g)
#         locality_init = loc_1step!(P_initial, C, N_side), loc_1step_0!(P_initial, C, N_side)
#         locality_final = loc_1step!(P, C, N_side), loc_1step_0!(P, C, N_side)
#         energy_final = energy(g, P, C)
#         SA_extremal = P_initial, energy_initial, nr_gen_con(P_initial, N_side), P, energy_final, nr_gen_con(P, N_side), N_T, en, locality_init, locality_final
#         push!(Data, SA_extremal)
#     end
#     Data
# end

### data collection: collects multiple runs of eval_sim_anneal!() in one object
function collect_data_SA_runs(N_runs, N_side, C, T, k_max)
    Data = []
    for i in 1:N_runs
        SA_extremal = eval_sim_anneal!(N_side, C, T, 0, k_max)
        push!(Data, SA_extremal)
    end
    Data
end

### safe data
using JLD
using Dates ##### ToDo for saving date and time in filename
#Dates.now(Dates.UTC)

### visualize data
function visualize_data(SData, rand_opt, Run_Nr) # for random grid: rand_opt=1 for optimized grid: rand_opt=4
    Data = SData["Data"]
    N_side = SData["N_side"]
    P = Data[Run_Nr][rand_opt]
    g = gen_square_grid(N_side)
    visualize_graph(g, P, N_side)
end

function flows_above_thres(T, P_SA, g) # gives number of flows above threshold T
    F = flow(g, P_SA)
    count(x -> x > T, abs.(F))
end

function collect_data(Data, col)
    Data_x = [ ]
    N_runs = length(Data)
    for i in 1:N_runs
        x = Data[i][col]
        push!(Data_x, x)
    end
    Data_x
end

function collect_data2(Data, col, subcol)
    Data_x = [ ]
    N_runs = length(Data)
    for i in 1:N_runs
        x = Data[i][col][subcol]
        push!(Data_x, x)
    end
    Data_x
end

function collect_av_e_std(SData)
    Data = SData["Data"]
    N_runs = SData["N_runs"]
    Data_mean = [ ]
    Data_std = [ ]
    N_steps = length(Data[1][8])
    for i in 1:N_steps
        x = mean(collect_data2(Data, 8, i))
        push!(Data_mean, x)
        y = 1 / sqrt(N_runs) * std(collect_data2(Data, 8, i); corrected=true)
        push!(Data_std, y)
    end
    Data_mean, Data_std
end

function collect_av_e_diff_std(Data)
    Data_mean = [ ]
    Data_std = [ ]
    N_runs = length(Data[1][6])
    for i in 2:N_runs
        x = mean(collect_data2(Data, 6, i - 1)) - mean(collect_data2(Data, 6, i))
        push!(Data_mean, x)
        y = sqrt((1 / sqrt(N_runs) * std(collect_data2(Data, 6, i - 1); corrected=true))^2 + (1 / sqrt(N_runs) * std(collect_data2(Data, 6, i); corrected=true))^2)
        push!(Data_std, y)
    end
    Data_mean, Data_std
end

function energy_from_data(Data, N_runs, C, N_side) # davor energy() abändern
    edge_energy = [ ]
    for i in 1:N_runs
        g = gen_square_grid(N_side)
        energy = energy(g, Data[i][3], C, N_side)
        push!(edge_energy, energy[2])
    end
    #append!(Data, edge_energy)
    #Data
    edge_energy
end

function energy_from_data2(Data, N_runs, C, N_side)
    g_av_energy = [ ]
    for i in 1:N_runs
        energy = Data[i][4][2]
        push!(g_av_energy, energy)
    end
    #append!(Data, edge_energy)
    #Data
    g_av_energy
end
