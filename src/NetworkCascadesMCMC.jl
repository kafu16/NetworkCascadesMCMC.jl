module NetworkCascadesMCMC

#### ToDo: export functions
#= dont export funcs that aren't made for interactive use and are used only
locally in a script =#

using LinearAlgebra

include("network_topologies.jl")
include("power_injections.jl")
include("flow_calculation.jl")
include("simulated_annealing.jl")
include("helpers.jl")
include("visualization.jl")
include("postprocessing.jl")
include("visualization.jl")

using Plots
### export of plots
using Compose
using Cairo
using Fontconfig
#
# # save to pdf
# draw(PDF("bla.pdf", 16cm, 16cm), gplot(g))
# # save to png
# draw(PNG("bla.png", 16cm, 16cm), gplot(g))
# # save to svg
# draw(SVG("bla_neu.svg", 16cm, 16cm), gplot(g))

# performance
using CPUTime

# Hier nur Hauptfunktion, die alles zusammensetzt.

# standard parameters
N_side = 4
C = 1
Steps_per_temp = 4 # number of steps before lowering temperature
k_max = 10
ann_shed = "0.999 ^ k"
N_runs = 1
T = 0.95
Run_Nr = 1

Data = collect_data_SA_runs_var_ann_shed(N_runs, N_side, C, T, temp_ex1, k_max)

end
