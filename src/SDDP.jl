#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module SDDP

import Reexport
Reexport.@reexport using JuMP

import JSON, Printf, Random, RecipesBase, TimerOutputs, Statistics
import MathOptFormat

export @stageobjective

# Modelling interface.
include("user_interface.jl")

# Default definitions for SDDP related modular utilities.
include("plugins/headers.jl")

# Printing utilities.
include("print.jl")

# The core SDDP code.
include("algorithm.jl")

# Specific plugins.
include("plugins/risk_measures.jl")
include("plugins/sampling_schemes.jl")
include("plugins/bellman_functions.jl")
include("plugins/stopping_rules.jl")

# Visualization related code.
include("visualization/publication_plot.jl")
include("visualization/spaghetti_plot.jl")

end
