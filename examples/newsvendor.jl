#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

srand(11111)

# Demand for newspapers
# There are two equally probable noises in each stage
#   Demand[stage, noise]
Demand = [
    10. 15.;
    12. 20.;
    8.  20.
]

# Markov state purchase prices
PurchasePrice = [5., 8.]

RetailPrice = 7.

# Transition matrix
Transition = Array{Float64, 2}[
    [1.0]',
    [0.6 0.4],
    [0.3 0.7; 0.3 0.7]
  ]

# Initialise SDDP Model
m = SDDPModel(
        sense             = :Min,
        stages            = 3,
        objective_bound   = -1000,
        markov_transition = Transition,
        solver            = ClpSolver(),
        cut_oracle        = DematosCutOracle()
                                                ) do sp, stage, markov_state

    # ====================
    #   State variable
    @state(sp, 0 <= stock <= 100, stock0==5)

    # ====================
    #   Other variables
    @variables(sp, begin
        buy  >= 0  # Quantity to buy
        sell >= 0  # Quantity to sell
    end)

    # ====================
    #   Noises
    @noises(sp, D=Demand[stage,:], begin
        sell <= D
        sell >= 0.5D
    end)

    # ====================
    #   Objective
    stageobjective!(sp, -sell * RetailPrice + buy * PurchasePrice[markov_state])

    # ====================
    #   Dynamics constraint
    @constraint(sp, stock == stock0 + buy - sell)

end

@time solvestatus = SDDP.solve(m,
    max_iterations = 20,
    cut_selection_frequency = 10,
    simulation     = MonteCarloSimulation(
                        frequency = 10,
                        min       = 10,
                        max       = 500,
                        step      = 10
                             ),
    bound_convergence = BoundConvergence(
                        iterations = 5,
                        atol       = 1e-3
                            )
)

@test isapprox(getbound(m), -97.9, atol=1e-3)
@test solvestatus == :bound_convergence
