module Kokako

using Reexport
@reexport using JuMP

using JSON, Logging

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

end
