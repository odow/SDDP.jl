#  Copyright 2017, Oscar Dowson and contributors
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

function getmodel(;oracle=DefaultCutOracle())
    # Initialise SDDP Model
    m = SDDPModel(
            sense             = :Max,
            stages            = 3,
            objective_bound   = 1000,
            markov_transition = Transition,
            cut_oracle        = oracle,
            risk_measure      = NestedAVaR(
                                 beta   = 0.6,
                                 lambda = 0.5
                                 ),
            solver          = ClpSolver()
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
        stageobjective!(sp, sell * RetailPrice - buy * PurchasePrice[markov_state])

        # ====================
        #   Dynamics constraint
        @constraint(sp, stock == stock0 + buy - sell)

    end
    m
end

# Build and solve a model
m1 = getmodel()
cuts_file_name = "newsvendor_cuts.csv"
solvestatus = SDDP.solve(m1,
    max_iterations = 30,
    cut_output_file = cuts_file_name)

@assert solvestatus == :max_iterations

srand(22222)
results = simulate(m1, 500)
@test isapprox(mean(r[:objective] for r in results), 97.92, atol=0.1)

# Build a completely new model and load old cuts
m2 = getmodel()
loadcuts!(m2, cuts_file_name)

# Simulating the new model gives the same results as before
srand(22222)
results = simulate(m2, 500)
@test isapprox(mean(r[:objective] for r in results), 97.92, atol=0.1)

# Build a completely new model with a different cut manager
m3 = getmodel(oracle=DematosCutOracle())
# Even dominated cuts are added and kept
loadcuts!(m3, cuts_file_name)

# Check again that we have the same/very similar results
srand(22222)
results = simulate(m3, 500)
@test isapprox(mean(r[:objective] for r in results), 97.92, atol=0.1)

rm(cuts_file_name)
