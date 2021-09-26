module NetworkCascadesMCMC

####################### import external packages ###############################
using LinearAlgebra

# performance
# sing CPUTime

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
export collect_data_SA_runs_var_ann_shed, collect_data_SA_runs



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
