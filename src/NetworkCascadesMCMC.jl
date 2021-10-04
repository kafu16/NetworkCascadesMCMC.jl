module NetworkCascadesMCMC

####################### import external packages ###############################
using LinearAlgebra

############################ source files ######################################
include("network_topologies.jl")
include("power_injections.jl")
include("flow_calculation.jl")
include("simulated_annealing.jl")
include("helpers.jl")
include("visualization.jl")
include("postprocessing.jl")

########################### export functions ###################################
# core code simulated annealing
export energy, sim_anneal, multiple_sim_anneal, parallel_multiple_sim_anneal
export flow, linefailure!, cascade!

# postprocessing
export postprocess_sim_anneal, postprocess_sim_anneal_high_gc_low_Gav
export write_out_params, write_out_postprocess, write_out_energies
export flows_above_thres, nr_gen_con, nr_gen_con_av, locality

# plotting
export plot_Gav_single_run, plot_Gav_av
export plot_histogram_all_runs, plot_histogram_high_gc_low_Gav

# visualization/plotting
export visualize_graph, visualize_data, visualize_graph_after_linefailure_cascade
export plot_Gav_av, plot_Gav_single_run, plot_histogram_random_vs_minimized_G_av,
plot_histogram_all_runs, plot_histogram_high_gc_low_Gav

# square grids
export gen_square_grid, gen_periodic_square_grid
export gen_stable_square_config, gen_multiple_stable_square_configs

# power injections
export gen_stable_square_config, gen_multiple_stable_square_configs

# annealing schedules
export temp_ex1_a, temp_ex1_b, temp_ex1_c, temp_ex2_a, temp_ex3_a, temp_ex3_b, temp_ex4_a, temp_ex5_a

# old functions
export collect_data_SA_runs


# # save to pdf
# draw(PDF("bla.pdf", 16cm, 16cm), gplot(g))
# # save to png
# draw(PNG("bla.png", 16cm, 16cm), gplot(g))
# # save to svg
# draw(SVG("bla.svg", 16cm, 16cm), gplot(g))

end


#= ToDo
 - [ ] define types for function arguments ([2021-08-20 Fr] done for simulated_annealing.jl)
 - [ ] Performance: code linefailure!() with sparse matrix
 - [ ] Performance: instead of calculating incidence matrix multiple time anew actualize entries of the matrix
 - [ ] if advantageous do package internal data management using DataFrames
 - [ ] visualization: heat map for flows
=#
