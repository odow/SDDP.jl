#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

# For repeatability
srand(11111)

include(joinpath(dirname(@__FILE__), "hydro_valley.jl"))

# deterministic
deterministic_model = hydrovalleymodel(hasmarkovprice=false,hasstagewiseinflows=false)
SDDP.solve(deterministic_model, max_iterations=10, cut_selection_frequency=1)
@test isapprox(getbound(deterministic_model), 835.0, atol=1e-3)

# stagewise inflows
stagewise_model = hydrovalleymodel(hasmarkovprice=false)
SDDP.solve(stagewise_model, max_iterations=20, print_level=0)
@test isapprox(getbound(stagewise_model), 838.33, atol=1e-2)

# markov prices
markov_model = hydrovalleymodel(hasstagewiseinflows=false)
SDDP.solve(markov_model, max_iterations=10, print_level=0)
@test isapprox(getbound(markov_model), 851.8, atol=1e-2)

# stagewise inflows and markov prices
markov_stagewise_model = hydrovalleymodel(hasstagewiseinflows=true, hasmarkovprice=true)
SDDP.solve(markov_stagewise_model, max_iterations=10, print_level=0)
@test isapprox(getbound(markov_stagewise_model), 855.0, atol=1e-3)

# risk averse stagewise inflows and markov prices
riskaverse_model = hydrovalleymodel(riskmeasure = NestedAVaR(lambda=0.5, beta=0.66))
SDDP.solve(riskaverse_model,
    max_iterations = 10,
    print_level = 0
)

@test isapprox(getbound(riskaverse_model), 828.157, atol=1e-3)

# stagewise inflows and markov prices
worst_case_model = hydrovalleymodel(
    riskmeasure = NestedAVaR(lambda=0.5, beta=0.0), sense=:Min)
SDDP.solve(worst_case_model,
    max_iterations = 10,
    simulation = MonteCarloSimulation(
        frequency = 2,
        min  = 20,
        step = 10,
        max  = 50
    )
)

@test isapprox(getbound(worst_case_model), -780.867, atol=1e-3)

# stagewise inflows and markov prices
cutselection_model = hydrovalleymodel(cutoracle = LevelOneCutOracle(), sense=:Max)
SDDP.solve(cutselection_model,
    max_iterations          = 10,
    print_level = 0,
    cut_selection_frequency = 2
)

@test length(SDDP.getsubproblem(cutselection_model, 2, 1).linconstr) < 10 + 11
@test isapprox(getbound(cutselection_model), 855.0, atol=1e-3)
