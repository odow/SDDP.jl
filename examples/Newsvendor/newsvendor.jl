#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
This problem is a the classic newsvendor problem. The agent has a stock of
newspapers to sell over a period of time, and needs to determine the quantity to
restock each day. However, demand is random and the price varies according to a
markov chain.
=#
using SDDP, JuMP, Clp, Base.Test

function newsvendormodel(;oracle=DefaultCutOracle(), riskmeasure=Expectation())
    # Demand for newspapers
    # There are two equally probable noises in each stage
    #   Demand[stage, noise]
    Demand = [
        10.0 15.0;
        12.0 20.0;
         8.0 20.0
    ]

    # Markov state purchase prices
    PurchasePrice = [5.0, 8.0]

    RetailPrice = 7.0

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
            cut_oracle        = oracle,
            risk_measure      = riskmeasure
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
        @rhsnoises(sp, D=Demand[stage,:], begin
            sell <= D
            sell >= 0.5D
        end)

        # ====================
        #   Objective
        @stageobjective(sp, -sell * RetailPrice + buy * PurchasePrice[markov_state])

        # ====================
        #   Dynamics constraint
        @constraint(sp, stock == stock0 + buy - sell)

    end
end

news1 = newsvendormodel()

@test SDDP.solve(news1,
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
) == :bound_convergence
@test isapprox(getbound(news1), -97.9, atol=1e-3)

news2 = newsvendormodel(riskmeasure=NestedAVaR(beta=0.6,lambda=0.5))
@test SDDP.solve(news2,
    max_iterations = 50,
    print_level = 0,
    simulation     = MonteCarloSimulation(
                        frequency = 10,
                        min       = 500,
                        max       = 500,
                        step      = 1
                             )
) == :max_iterations
@test isapprox(getbound(news2), -93.267, atol=1e-3)
results2 = simulate(news2, 500)
@test length(results2) == 500
@test isapprox(mean(r[:objective] for r in results2), -97, atol=2)


historical_results = simulate(news2, [:buy, :sell];
    markovstates = [1, 2, 1],
    noises    = [1, 1, 1]
)
@test isapprox(historical_results[:objective], -85)

historical_results2 = simulate(news2, [:buy, :sell];
    markovstates = [1, 2, 1],
    noises    = [2, 2, 2]
)
@test isapprox(historical_results2[:objective], -119)


# Build and solve a model
news3 = newsvendormodel(riskmeasure=NestedAVaR(beta=0.6,lambda=0.5))
cuts_file_name = "newsvendor_cuts.csv"
@test SDDP.solve(news3,
    max_iterations = 30,
    print_level = 0,
    cut_output_file = cuts_file_name) == :max_iterations

srand(22222)
results3 = simulate(news3, 500)
@test isapprox(mean(r[:objective] for r in results3), -97.9, atol=0.1)

# Build a completely new model and load old cuts
news4 = newsvendormodel(riskmeasure=NestedAVaR(beta=0.6,lambda=0.5))
loadcuts!(news4, cuts_file_name)

# Simulating the new model gives the same results as before
srand(22222)
results4 = simulate(news4, 500)
@test isapprox(mean(r[:objective] for r in results4), -97.9, atol=0.1)

# Build a completely new model with a different cut manager
news5 = newsvendormodel(riskmeasure=NestedAVaR(beta=0.6,lambda=0.5), oracle=DematosCutOracle())
# Even dominated cuts are added and kept
loadcuts!(news5, cuts_file_name)

# Check again that we have the same/very similar results
srand(22222)
results5 = simulate(news5, 500)
@test isapprox(mean(r[:objective] for r in results5), -97.9, atol=0.1)

rm(cuts_file_name)
