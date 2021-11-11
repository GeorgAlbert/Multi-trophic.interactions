############################################################################
# Framework for creating and simulating the dynamics of food webs used in  #
# Albert et al. (202x): The hidden role of multi-trophic interactions in   #
# driving diversity-productivity relationships. Ecology Letters            #
############################################################################

## General remarks
# simulation depends on specific intrinsic data-structure, i.e. 1) order of parameter-vector is always sorted the same, 2) if animals, producers, and/or resources have shared parameters, values are ordered as [animals..., producers..., resources...], and 3) whenever there is a species x species matrix (e.g. for feeding rates), rows are resource species and columns are consumer species
# for feedback or questions, please contact: georg.albert.ecol@gmail.com
# for the most up to date version, visit: https://github.com/GeorgAlbert/Multi-trophic.interactions


## load required packages and functions
using DifferentialEquations
using Distributions

include("functions.jl")


## example simulation

# 1. generate initial set of starting parameters and densities
# - 30 animal and 16 producer species:
start = startparam(30, 16)
# - with scenario of RUD gradient (use index to specify which scenario, where 1 represents RUD = 1 and the max value represents RUD = 0)
start = startparam(30, 16, nC = 16, RUDscenario = RUD(16, 16)[1])
# - with random RUD scenario
start = startparam(30, 16, nC = 16, RUDscenario = RUD(16, 16, random = true))


# 2. modify staring parameters and densities
# - replace RUD scenario while keeping everything else constant
start = modify_startparam(start, RUD = RUD(16, 16)[16])
# - select producer species to keep and remove the rest
start = modify_startparam(start, [1, 2])
# - replace animal community, e.g. to reduce their diversity
start = modify_startparam(start, nA_new = 10)


# 3. run simulation
# - create ODE problem
tspan = (0.0, 150000.0)
prob = ODEProblem(foodwebsim, start[2], tspan, start[1])
# - solve ODE problem
sol = solve(prob, AutoVern7(Rodas4()), reltol = 1e-8, abstol = 1e-8, maxiters = 15000)
