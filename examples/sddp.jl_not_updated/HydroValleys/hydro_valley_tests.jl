#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Test, Random

# For repeatability
Random.seed!(11111)

include(joinpath(dirname(@__FILE__), "hydro_valley.jl"))

# deterministic
deterministic_model = hydrovalleymodel(hasmarkovprice=false,hasstagewiseinflows=false)
SDDP.train(deterministic_model, iteration_limit=10, cut_deletion_minimum=1, print_level=0)
@test SDDP.calculate_bound(deterministic_model) ≈ 835.0 atol = 1e-3

# stagewise inflows
stagewise_model = hydrovalleymodel(hasmarkovprice=false)
SDDP.train(stagewise_model, iteration_limit=20, print_level=0)
@test SDDP.calculate_bound(stagewise_model) ≈ 838.33 atol = 1e-2

# markov prices
markov_model = hydrovalleymodel(hasstagewiseinflows=false)
SDDP.train(markov_model, iteration_limit=10, print_level=0)
@test SDDP.calculate_bound(markov_model) ≈ 851.8 atol = 1e-2

# stagewise inflows and markov prices
markov_stagewise_model = hydrovalleymodel(hasstagewiseinflows=true, hasmarkovprice=true)
SDDP.train(markov_stagewise_model, iteration_limit=10, print_level=0)
@test SDDP.calculate_bound(markov_stagewise_model) ≈ 855.0 atol = 1e-3

# risk averse stagewise inflows and markov prices
riskaverse_model = hydrovalleymodel()
SDDP.train(riskaverse_model, risk_measure = SDDP.EAVaR(lambda=0.5, beta=0.66), iteration_limit = 10, print_level = 0)
@test SDDP.calculate_bound(riskaverse_model) ≈ 828.157 atol = 1e-3

# stagewise inflows and markov prices
worst_case_model = hydrovalleymodel(sense=:Min)
SDDP.train(
    worst_case_model
    , risk_measure = SDDP.EAVaR(lambda=0.5, beta=0.0)
    , iteration_limit = 10
    , print_level=0
)
@test SDDP.calculate_bound(worst_case_model) ≈ -780.867 atol = 1e-3

# stagewise inflows and markov prices
cutselection_model = hydrovalleymodel()
SDDP.train(
    cutselection_model
    , iteration_limit = 10
    , print_level = 0
    , cut_deletion_minimum = 2
)
# @test length(SDDP.getsubproblem(cutselection_model, 2, 1).linconstr) < 10 + 11
@test SDDP.calculate_bound(cutselection_model) ≈ 855.0 atol = 1e-3

# Distributionally robust Optimization
dro_model = hydrovalleymodel(hasmarkovprice = false)
SDDP.train(dro_model, risk_measure = SDDP.ModifiedChiSquared(sqrt(2/3) - 1e-6), iteration_limit = 10, print_level = 0)
@test SDDP.calculate_bound(dro_model) ≈ 835.0 atol = 1e-3

dro_model = hydrovalleymodel(hasmarkovprice = false)
SDDP.train(dro_model, risk_measure = SDDP.ModifiedChiSquared(1/6), iteration_limit = 20, print_level = 0)
@test SDDP.calculate_bound(dro_model) ≈ 836.695 atol = 1e-3
# (Note) radius ≈ sqrt(2/3), will set all noise probabilities to zero except the worst case noise
# (Why?):
# The distance from the uniform distribution (the assumed "true" distribution)
# to a corner of a unit simplex is sqrt(S-1)/sqrt(S) if we have S scenarios. The corner
# of a unit simplex is just a unit vector, i.e.: [0 ... 0 1 0 ... 0]. With this probability
# vector, only one noise has a non-zero probablity.
# In the worst case rhsnoise (0 inflows) the profit is:
#  Reservoir1: 70 * $3 + 70 * $2 + 65 * $1 +
#  Reservoir2: 70 * $3 + 70 * $2 + 70 * $1
#  = $835
