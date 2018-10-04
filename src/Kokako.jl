module Kokako

import Reexport
Reexport.@reexport using JuMP

import JSON, Printf, RecipesBase, Statistics

export @stageobjective

# Modelling interface.
include("user_interface.jl")

# SDDP related modular utilities.
include("risk_measures.jl")
include("sampling_schemes.jl")
include("bellman_functions.jl")
include("stopping_rules.jl")

# Printing utilities.
include("print.jl")

# The core SDDP code.
include("sddp.jl")

include("policy_visualization.jl")

end
