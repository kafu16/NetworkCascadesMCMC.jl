module NetworkCascadesMCMC

####################### import external packages ###############################
using LinearAlgebra
using Dates

# performance
using CPUTime

############################ source files ######################################
include("network_topologies.jl")
include("power_injections.jl")
include("flow_calculation.jl")
include("simulated_annealing.jl")
include("helpers.jl")
include("visualization.jl")
include("postprocessing.jl")

########################### export functions ###################################
#### ToDo: export functions
#= dont export funcs that aren't made for interactive use and are used only
locally in a script =#


# core code simulated annealing
export energy, sim_anneal, multiple_sim_anneal

# postprocessing
export postprocess_sim_anneal, postprocess_sim_anneal_high_gc_low_Gav
export write_out_params, write_out_postprocess

# plotting
export plot_Gav_single_run, plot_Gav_av
export plot_histogram_all_runs, plot_histogram_high_gc_low_Gav

# visualization
export visualize_graph, visualize_data

# square grids
export gen_square_grid
export gen_stable_square_config, gen_multiple_stable_square_configs

export temp_ex5

# old
export collect_data_SA_runs_var_ann_shed, collect_data_SA_runs,

# Hier nur Hauptfunktion, die alles zusammensetzt.


################################################################################
########################## DELETE AT END OF SESSION ############################
################################################################################

############################ PARAMETERS ########################################
N_side = 4
C = 1.
steps_per_temp = 4 # number of steps before lowering temperature
k_max = 100
#k_max = 2291
annealing_schedule = temp_ex5
ann_sched = "0.99^k" # for storing annealing schedule in .jld
N_runs = 4
T = 0.95
Run_Nr = 1
linenumber = 1
gen_con = 11
G_av_final = 10
################################################################################

function main()
    repo_directory = pwd()
    t=now()
    datetime = Dates.format(t, "yyyymmdd_HHMMSS.s") # https://riptutorial.com/julia-lang/example/20476/current-time
    folder = string("data/",datetime,"_N_runs",string(N_runs),"_k_max",string(k_max),"_ann_sched",string(ann_sched))
    directory = string(repo_directory,"/data/",datetime,"_N_runs",string(N_runs),"_k_max",string(k_max),"_ann_sched",string(ann_sched))
    mkpath(folder)
    cd(directory)

    g = gen_square_grid(N_side)
    P_inits = gen_multiple_stable_square_configs(g, N_side, C, N_runs)
    simulation_data = string(directory,"/simulation_data.jld")
    postprocess_data = string(directory,"/postprocess_data.jld")
    multiple_sim_anneal(simulation_data, g, P_inits, C, annealing_schedule, steps_per_temp, k_max)
    postprocess_sim_anneal(simulation_data, postprocess_data, T)
    Data_loaded = JLD.load(postprocess_data)

    filename = "params_postprocess.txt"
    write_out_params(Data_loaded, filename)
    write_out_postprocess(Data_loaded, filename)

    plot_Gav_av(Data_loaded)
    plot_histogram_all_runs(Data_loaded)
    Run_Nr = 1; plot_Gav_single_run(Data_loaded, Run_Nr)

    postprocess_high_gc_low_Gav_data = string(directory,"/postprocess_high_gc_low_Gav_data.jld")
    postprocess_sim_anneal_high_gc_low_Gav(simulation_data, postprocess_high_gc_low_Gav_data, T, gen_con, G_av_final)
    Data_loaded_high_gc_low_Gav = JLD.load(postprocess_high_gc_low_Gav_data)
    filename = "params_postprocess_high_gc_low_Gav.txt"
    write_out_params(Data_loaded_high_gc_low_Gav, filename)
    write_out_postprocess(Data_loaded_high_gc_low_Gav, filename)
    plot_histogram_high_gc_low_Gav(Data_loaded, Data_loaded_high_gc_low_Gav)

    cd(repo_directory)
end

main()

# collect_data_SA_runs_var_ann_shed()
Data = collect_data_SA_runs(N_runs, N_side, C, T, annealing_schedule, steps_per_temp, k_max)

g = gen_square_grid(N_side)
F_coll = flow(g, Data[1][4])
maximum(abs.(F_coll))

# sim_anneal()
g = gen_square_grid(N_side)
P_init = gen_stable_square_config(g, N_side, C)
P_sim, en = sim_anneal(g, P_init, C, temp_ex5, steps_per_temp, k_max)
energy(g, P_sim, C)

F_sim = flow(g, P_sim)
maximum(abs.(F_sim))

visualize_graph(g,P,N_side)
################################################################################
################################################################################
################################################################################


# # save to pdf
# draw(PDF("bla.pdf", 16cm, 16cm), gplot(g))
# # save to png
# draw(PNG("bla.png", 16cm, 16cm), gplot(g))
# # save to svg
# draw(SVG("bla.svg", 16cm, 16cm), gplot(g))

end


#= ToDo
 - [ ] define types for function arguments ([2021-08-20 Fr] done for simulated_annealing.jl)
=#
