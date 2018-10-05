module Kokako

import Reexport
Reexport.@reexport using JuMP

import JSON, Printf, RecipesBase, TimerOutputs, Statistics

export @stageobjective

# Modelling interface.
include("user_interface.jl")

# SDDP related modular utilities.
include("plugins/risk_measures.jl")
include("plugins/sampling_schemes.jl")
include("plugins/bellman_functions.jl")
include("plugins/stopping_rules.jl")

# Printing utilities.
include("print.jl")

# The core SDDP code.
include("sddp.jl")

include("policy_visualization.jl")

end
