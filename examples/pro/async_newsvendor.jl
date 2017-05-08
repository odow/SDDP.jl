#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

srand(11111)

# Demand for newspapers
# There are two equally probable scenarios in each stage
#   Demand[stage, scenario]
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
        sense             = :Max,
        stages            = 3,
        objective_bound   = 1000,
        markov_transition = Transition,
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
    #   Scenarios
    @scenarios(sp, D=Demand[stage,:], begin
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

@time solvestatus = SDDP.solve(m,
    solve_type     = Asyncronous(),
    max_iterations = 50,
    simulation     = MonteCarloSimulation(
                        frequency = 10,
                        min       = 5,
                        max       = 50,
                        step      = 5
                             )
)

@test isapprox(getbound(m), 93.267, atol=1e-3)
@test solvestatus == :max_iterations

@everywhere begin
    if myid() != 1
        mm = SDDP.m
        sp = SDDP.getsubproblem(mm, 1, 1)
        @show length(sp.linconstr)
        sp2 = SDDP.getsubproblem(mm, 2, 1)
        @show length(sp2.linconstr)
    end
end
sp = SDDP.getsubproblem(m, 1, 1)
@show length(sp.linconstr)
sp2 = SDDP.getsubproblem(m, 2, 1)
@show length(sp2.linconstr)
