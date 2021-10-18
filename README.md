# NetworkCascadesMCMC

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kafu16.github.io/NetworkCascadesMCMC.jl/)

## Installation
 - Install NetworkCascadesMCMC.jl: `pkg> add https://github.com/kafu16/NetworkCascadesMCMC.jl.git`
 - Get recent version of NetworkCascadesMCMC.jl: `pkg> update NetworkCascadesMCMC`

## Example Usage
This is a simple example script
```julia
# creating data folder
print("create folder and paths: ")
A = Dates.now()

repo_directory = pwd()
t=now()
datetime = Dates.format(t, "yyyymmdd_HHMMSS.s") # https://riptutorial.com/julia-lang/example/20476/current-time
folder = string("data/",datetime,"_N_runs",string(N_runs),"_k_max",string(k_max),"_ann_sched",string(annealing_schedule_name))
directory = string(repo_directory,"/data/",datetime,"_N_runs",string(N_runs),"_k_max",string(k_max),"_ann_sched",string(annealing_schedule_name))
mkpath(folder)
cd(directory)

B = Dates.now()
elapsed_time(A, B)


# core simulation
print("\ncore simulation: ")
A = Dates.now()
```
