#= this file contains helper functions that do not fit in any other file of
NetworkCascadesMCMC including plotting functions=#

using Plots
using StatsPlots, StatsBase
using LaTeXStrings

# export of plots
using Compose
using Cairo
using Fontconfig

#safe data
using JLD

# data collection sqaure grids
function collect_data_SA_runs(N_runs, N_side, C, T, annealing_schedule, steps_per_temp, k_max)
    Data = []
    for i in 1:N_runs
        g = gen_square_grid(N_side)
        P_init = gen_stable_square_config(g, N_side, C) # to avoid iteration steps it is important to start with a stable configurations see comment at stable_swapped_config!()
        P, en = sim_anneal(g, P_init, C, annealing_schedule, steps_per_temp, k_max)

        # calculating observables
        energy_initial = energy(g, P_init, C)
        N_T = flows_above_thres(T, P_init, g), flows_above_thres(T, P, g)
        locality_init = loc_1step(g, P_init, C), loc_1step_0(g, P_init, C)
        locality_final = loc_1step(g, P, C), loc_1step_0(g, P, C)
        energy_final = energy(g, P, C)
        SA_extremal = P_init, energy_initial, nr_gen_con(g, P_init), P, energy_final, nr_gen_con(g, P), N_T, en, locality_init, locality_final
        push!(Data, SA_extremal)
    end
    Data
end

""" Writes parameters of either the core-simulation-.jld-file  or the
    postprocess-.jld-file to a .txt-file.
"""
function write_out_params(Data_loaded, filename)
    open(filename, "a") do io
        write(io, "PARAMETERS\n")
        write(io, "Grid: "); write(io, string(Data_loaded["Grid"])); write(io, "\n")
        write(io, "N_vertices: "); write(io, string(Data_loaded["N_vertices"])); write(io, "\n")
        write(io, "Annealing schedule: "); write(io, Data_loaded["annealing_schedule"]); write(io, "\n")
        write(io, "Steps per temp: "); write(io, string(Data_loaded["steps_per_temp"])); write(io, "\n")
        write(io, "k_max: "); write(io, string(Data_loaded["k_max"])); write(io, "\n")
        write(io, "C: "); write(io, string(Data_loaded["C"])); write(io, "\n")
        write(io, "N_runs: "); write(io, string(Data_loaded["N_runs"])); write(io, "\n")
        write(io, "\n")
    end
end

""" Writes postprocessing data of the postprocess-.jld-file to a .txt-file.
"""
function write_out_postprocess(Data_loaded, filename)
    gen_gen_av_init, gen_gen_std_init, con_con_av_init, con_con_std_init, gen_con_av_init, gen_con_std_init = nr_gen_con_av(Data_loaded, "nr_gen_con_init")
    gen_gen_av_final, gen_gen_std_final, con_con_av_final, con_con_std_final, gen_con_av_final, gen_con_std_final = nr_gen_con_av(Data_loaded, "nr_gen_con_final")
    loc = locality(Data_loaded)

    open(filename, "a") do io
        write(io, "POSTPROCESSING DATA\n\n")

        write(io, "Gav_av for random and optimized configurations:\n")
        write(io, "Gav_av_init: "); write(io, string(Gav_av_STD_Gav(Data_loaded)[1]))
        write(io, " +/- "); write(io, string(Gav_av_STD_Gav(Data_loaded)[2])); write(io, "\n")
        write(io, "Gav_av_final: "); write(io, string(Gav_av_STD_Gav(Data_loaded)[3]))
        write(io, " +/- "); write(io, string(Gav_av_STD_Gav(Data_loaded)[4])); write(io, "\n")
        write(io, "\n\n")

        write(io, "Number of edges between generators and consumers for random and optimized configurations:\n")
        write(io, "Random configurations\n")
        write(io, "gen-gen "); write(io, string(gen_gen_av_init)); write(io, " +/- "); write(io, string(gen_gen_std_init))
        write(io, "\ncon-con "); write(io, string(con_con_av_init)); write(io, " +/- "); write(io, string(con_con_std_init))
        write(io, "\ngen-con "); write(io, string(gen_con_av_init)); write(io, " +/- "); write(io, string(gen_con_std_init))
        write(io, "\n\n")
        write(io, "Configurations with optimized G_av\n")
        write(io, "gen-gen "); write(io, string(gen_gen_av_final)); write(io, " +/- "); write(io, string(gen_gen_std_final))
        write(io, "\ncon-con "); write(io, string(con_con_av_final)); write(io, " +/- "); write(io, string(con_con_std_final))
        write(io, "\ngen-con "); write(io, string(gen_con_av_final)); write(io, " +/- "); write(io, string(gen_con_std_final))
        write(io, "\n\n")

        write(io, "Locality:\n")
        write(io, "Random configurations\n")
        write(io, "<D_1> "); write(io, string(loc[1][5])); write(io, " +/- "); write(io, string(loc[1][6]))
        write(io, " <D_0> "); write(io, string(loc[2][5])); write(io, " +/- "); write(io, string(loc[1][6]))

        write(io, "\n\n")
        write(io, "Configurations with optimized G_av\n")
        write(io, "<D_1> "); write(io, string(loc[1][7])); write(io, " +/- "); write(io, string(loc[1][8]))
        write(io, " <D_0> "); write(io, string(loc[2][7])); write(io, " +/- "); write(io, string(loc[2][8]))
        write(io, "\n\n")
    end
end


function flows_above_thres(T, P, g) # gives number of flows above threshold T
    F = flow(g, P)
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

#### TODo is this function still used?
function collect_data2(Data, col, subcol)
    Data_x = [ ]
    N_runs = length(Data)
    for i in 1:N_runs
        x = Data[i][col][subcol]
        push!(Data_x, x)
    end
    Data_x
end

""" Calculates energy and standard deviation for every iteration step k averaged
    over N_runs.
"""
function collect_av_e_std(Data_loaded)
    N_runs = Data_loaded["N_runs"]
    k_max = Data_loaded["k_max"]
    energies = Data_loaded["energies"]
    Data_mean = [ ]
    Data_std = [ ]

    for i in 1:k_max
        # calculate mean for each iteration step k
        Data_x = [ ] # energy values for each run at iteration step k
        for j in 1:N_runs
            x = energies[j][i]
            push!(Data_x, x)
        end
        x = mean(Data_x)
        push!(Data_mean, x)

        y = 1 / sqrt(N_runs) * std(Data_x; corrected=true)
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

function postprocess_sim_anneal_high_gc_low_Gav(filepath_in, filepath_out, T, gen_con, G_av_final)
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

    energies_high_gc_low_Gav = []
    P_inits_high_gc_low_Gav = []
    P_finals_high_gc_low_Gav = []

    energy_init = []; energy_final = []
    N_T_init = []; N_T_final = []
    loc_1step_init = []; loc_1step_final = []
    loc_1step_0_init = []; loc_1step_0_final = []
    gen_gen_init = []; con_con_init = []; gen_con_init = []
    gen_gen_final = []; con_con_final = []; gen_con_final = []


    for i in 1:N_runs
        if nr_gen_con(g,P_finals[i])[3] > gen_con && energies[i][k_max] < G_av_final

            push!(energies_high_gc_low_Gav, Data_loaded["energies"][i])
            push!(P_inits_high_gc_low_Gav, Data_loaded["P_inits"][i])
            push!(P_finals_high_gc_low_Gav, Data_loaded["P_finals"][i])

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
    end

    N_runs_high_gc_low_Gav = length(energy_final)
    nr_gen_con_init = Nr_gen_con(gen_gen_init, con_con_init, gen_con_init)
    nr_gen_con_final = Nr_gen_con(gen_gen_final, con_con_final, gen_con_final)
    locality = Locality(loc_1step_init,loc_1step_final,loc_1step_0_init,loc_1step_0_final)

    JLD.save(filepath_out, "energies",energies_high_gc_low_Gav, "P_inits",P_inits_high_gc_low_Gav, "P_finals",P_finals_high_gc_low_Gav, "N_vertices",N_vertices, "Grid",g,
        "annealing_schedule",ann_sched, "steps_per_temp",steps_per_temp, "C",C , "k_max",k_max, "N_runs",N_runs_high_gc_low_Gav,
        "energy_init",energy_init, "energy_final",energy_final, "N_T_init",N_T_init, "N_T_final",N_T_final,
        "locality",locality, "nr_gen_con_init",nr_gen_con_init, "nr_gen_con_final",nr_gen_con_final)
end

################################################################################
################################## PLOTTING ####################################
################################################################################

function visualize_data(Data_loaded, P_rand_opt, Run_Nr) # for random grid: rand_opt=1 for optimized grid: rand_opt=4
    P = Data_loaded[P_rand_opt][Run_Nr]
    g = Data_loaded["Grid"]
    visualize_graph(g, P)
end

""" Plots energy and standard deviation for every iteration step k averaged over
    N_runs.
"""
function plot_Gav_av(Data_loaded)
    Data_av_e_std = collect_av_e_std(Data_loaded)
    Plots.plot(Data_av_e_std[1], title = L"\overline{G_{av}} \textrm{ depending on iteration step } k",
            label = L"\overline{G_{av}}",
            xaxis = L"\textrm{Iteration step } k", yaxis = L"\overline{G_{av}}",
            xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, titlefontsize = 10,
            ribbon=Data_av_e_std[2], fillalpha=.2,
            framestyle = :box, grid = true)
    #Plots.savefig("decreasing_Gav_av.pdf")
    Plots.savefig("decreasing_Gav_av")
end

""" Plots G_av for single run for each iteration step k.
"""
function plot_Gav_single_run(Data_loaded, Run_Nr)
    en = Data_loaded["energies"][Run_Nr]
    Plots.plot(en, title = L"\textrm{Decreasing } G_{av} \textrm{ for one sample grid}",
            label = L"G_{av}",
            xaxis = L"\textrm{Iteration step } k", yaxis = L"G_{av}",
            xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, titlefontsize = 10,
            framestyle = :box, grid = true)
    #Plots.savefig("decreasing_Gav_sample.pdf")
    Plots.savefig("decreasing_Gav_sample")
end

""" Calculates histogram for flows and calculates average and standard deviation
    for single bins, so the average and standard deviation over all first bins and
    separately for all second bins  etc. is calculated.
"""
function plot_histogram_random_vs_minimized_G_av(filename,random,minimized,nbins,xlabelstr,ylabelstr,labelstr,titlestr)
    weights1 = convert(Vector{Float64}, random[1])
    weights2 = convert(Vector{Float64}, minimized[1])
    errors1 = convert(Vector{Float64}, random[2])
    errors2 = convert(Vector{Float64}, minimized[2])

    # input in one vector where first half of vector refers to first group and second half to second group
    weights = append!(weights1, weights2)
    sx = repeat(["B", "A"], inner = nbins - 1)
    errors = append!(errors1, errors2)

    groupedbar(weights, yerr = errors, group = sx,
            xlabel = xlabelstr,
            xguidefontsize = 10, yguidefontsize = 10, legendfontsize = 8, titlefontsize = 10,
            ylabel = ylabelstr,
            label = labelstr,
            title = L"\textrm{Normalized histogram of flow distribution}",
            bar_width = 0.67,
            lw = 0.0, markerstrokewidth = 0.7, markerstrokecolor = :black,
            c = [:deepskyblue :orange], #https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
            framestyle = :box, grid = true, xticks = 1:1:nbins-1)
    Plots.savefig(filename)
end

function plot_histogram_all_runs(Data_loaded)
    nbins = 11
    random = weights_mean_err(Data_loaded, "P_inits", nbins)
    minimized = weights_mean_err(Data_loaded, "P_finals", nbins)
    # filename = "flows_rand_min_Gav.pdf"
    filename = "flows_rand_min_Gav"
    xlabelstr = L"\textrm{10 bins of width 0.10 from bin 1 = } [0.00,0.10) \textrm{ to bin 10 = } [0.90,1.00)"
    ylabelstr = L"\textrm{Normalized probability}"
    labelstr = [L"\textrm{Grids with low } G_{av}" L"\textrm{Random grids}"]
    titlestr = L"\textrm{Normalized histogram of flow distribution}"
    plot_histogram_random_vs_minimized_G_av(filename,random,minimized,nbins,xlabelstr,ylabelstr,labelstr,titlestr)
end

""" Compares the flows of ALL configurations to those configurations with high_gc_low_Gav.
"""
function plot_histogram_high_gc_low_Gav(Data_loaded, Data_loaded__high_gc_low_Gav)
    nbins = 21
    random = weights_mean_err(Data_loaded, "P_inits", nbins)
    minimized = weights_mean_err(Data_loaded__high_gc_low_Gav, "P_finals", nbins)
    # filename = "flows_rand_high_gc_low_Gav.pdf"
    filename = "flows_rand_high_gc_low_Gav"
    xlabelstr = L"\textrm{20 bins of width 0.05 from bin 1: } [0.00,0.05) \textrm{ to bin 20: } [0.95,1.00)"
    ylabelstr = L"\textrm{Normalized probability}"
    labelstr = [L"\textrm{Grids with high gc and low } G_{av}" L"\textrm{Random grids}"]
    titlestr = L"\textrm{Distribution: Mean of 1000 random and 17 grids with high gc and low } G_{av}"
    plot_histogram_random_vs_minimized_G_av(filename,random,minimized,nbins,xlabelstr,ylabelstr,labelstr,titlestr)
end

# histograms
function flow_single_sample(Data_loaded, i, P_rand_opt) # i: sample number
    g = Data_loaded["Grid"]
    P = Data_loaded[P_rand_opt][i]
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


""" Calculates histogram for flows and calculates average and standard deviation
    for single bins, so the average and standard deviation over all first bins and
    separately for all second bins  etc. is calculated.
"""
function weights_mean_err(Data_loaded, Config, nbins) # j = 1: random grids, j = 4: optimised grids
    weights = []
    N_bars = nbins - 1
    for i in 1:N_bars
        push!(weights, [])
    end

    N_runs = Data_loaded["N_runs"]
    bins = range(0.0, 1.0, length = nbins) # divides interval 0.0 and 1.0 in nbins sections
    for i in 1:N_runs
        flows = flow_single_sample(Data_loaded, i, Config)
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

################################################################################
######################### DEPRECATED FUNCTIONS##################################
################################################################################

# # data collection
# function collect_data_SA_runs_var_ann_shed(N_runs, N_side, C, T, annealing_schedule, k_max)
#     Data = []
#     for i in 1:N_runs
#         g = gen_square_grid(N_side)
#         P = gen_stable__square_config(g, N_side, C) # to avoid iteration steps it is important to start with a stable configurations see comment at stable_swapped_config!()
#         P_initial = copy(P)
#
#         en = [ ]
#         for k in 0:k_max - 1
#             Temp = annealing_schedule(k)
#             P_old = copy(P) # for calculating the energy of "old" configuration
#             P = stable_swapped_config!(g, P, C)
#             energy_old = energy(g, P_old, C)[2] # by [2] only the second value of tuple is returned (G_av)
#             energy_new = energy(g, P, C)[2] # by [2] only the second value of tuple is returned (G_av)
#             ΔE = energy_new - energy_old
#             #### performance: let energy() calculate G_av only
#
#             if ΔE <= 0 # man könnte conditional auch umdrehen: if (ΔE <= 0 AND probability(ΔE, T) < rand())
#                                                                      # P = P_old
#                 P
#             elseif probability(ΔE, Temp) > rand() # rand() gives random number element of [0,1]
#                 P
#             else
#                 P = P_old
#             end
#             g = gen_square_grid(N_side)
#             push!(en, energy_old)
#         end
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
