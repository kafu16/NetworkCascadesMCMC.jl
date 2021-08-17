module NetworkCascadesMCMC

#### ToDo: export functions
#= dont export funcs that aren't made for interactive use and are used only
locally in a script =#
export collect_data_SA_runs_var_ann_shed, collect_data_SA_runs

using LinearAlgebra

include("network_topologies.jl")
include("power_injections.jl")
include("flow_calculation.jl")
include("simulated_annealing.jl")
include("helpers.jl")
include("visualization.jl")
include("postprocessing.jl")

using Plots
### export of plots
using Compose
using Cairo
using Fontconfig

# # save to pdf
# draw(PDF("bla.pdf", 16cm, 16cm), gplot(g))
# # save to png
# draw(PNG("bla.png", 16cm, 16cm), gplot(g))
# # save to svg
# draw(SVG("bla.svg", 16cm, 16cm), gplot(g))

# performance
using CPUTime

# Hier nur Hauptfunktion, die alles zusammensetzt.


###############################################################################
###### DELETE AT END OF SESSION ###############################################
###############################################################################
N_side = 4
C = 1
Steps_per_temp = 1 # number of steps before lowering temperature
#k_max = 10
k_max = 2291
ann_shed = "0.99 ^ k"
N_runs = 1
T = 0.95
Run_Nr = 1

# collect_data_SA_runs_var_ann_shed()
Data = collect_data_SA_runs(N_runs, N_side, C, T, temp_ex5, k_max)
g = gen_square_grid(N_side)
F_coll = flow(g, Data[1][4])
maximum(abs.(F_coll))

# sim_anneal()
g = gen_square_grid(N_side)
P = gen_stable_config(g, N_side, C)
P_sim, en = sim_anneal(g, P, C, temp_ex5, k_max)
energy(g, P_sim, C)

F_sim = flow(g, P_sim)
maximum(abs.(F_sim))

###############################################################################
###############################################################################
###############################################################################


end
