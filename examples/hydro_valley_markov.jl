#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

# For repeatability
srand(11111)

immutable TurbineC
    flowknots::Vector{Float64}
    powerknots::Vector{Float64}
end

immutable ReservoirC
    min::Float64
    max::Float64
    initial::Float64
    turbine::TurbineC
    spill_cost::Float64
    inflows::Vector{Float64}
end

valley_chain = [
    ReservoirC(0, 200, 200, TurbineC([50, 60, 70], [55, 65, 70]), 1000, [0, 20, 50]),
    ReservoirC(0, 200, 200, TurbineC([50, 60, 70], [55, 65, 70]), 1000, [0, 0,  20])
]
turbine(i) = valley_chain[i].turbine

# Prices[stage, markov state]
prices = [
    1 2;
    2 1;
    3 4
]
# prices = [
#     1 1;
#     2 2;
#     3 3
# ]

# Transition matrix
transition = Array{Float64, 2}[
    [ 1.0 ]',
    [ 0.6 0.4 ],
    [ 0.6 0.4; 0.3 0.7 ]
]

N = length(valley_chain)

# Initialise SDDP Model
m = SDDPModel(
            sense           = :Max,
            stages          = 3,
            objective_bound = 1e6,
            markov_transition = transition,
            risk_measure    = Expectation(),
            cut_oracle      = DefaultCutOracle(),
            solver          = ClpSolver()
                                    ) do sp, stage, markov_state

    # ------------------------------------------------------------------
    #   SDDP State Variables
    # Level of upper reservoir
    @state(sp, valley_chain[r].min <= reservoir[r=1:N] <= valley_chain[r].max, reservoir0==valley_chain[r].initial)

    # ------------------------------------------------------------------
    #   Additional variables
    @variables(sp, begin
        outflow[r=1:N]      >= 0
        spill[r=1:N]        >= 0
        inflow[r=1:N]       >= 0
        generation_quantity >= 0 # Total quantity of water
        # Proportion of levels to dispatch on
        0 <= dispatch[r=1:N, level=1:length(turbine(r).flowknots)] <= 1
    end)

    # ------------------------------------------------------------------
    # Constraints
    @constraints(sp, begin
        # flow from upper reservoir
        reservoir[1] == reservoir0[1] + inflow[1] - outflow[1] - spill[1]
        # other flows
        flow[i=2:N], reservoir[i] == reservoir0[i] + inflow[i] - outflow[i] - spill[i] + outflow[i-1] + spill[i-1]

        # Total quantity generated
        generation_quantity == sum(turbine(r).powerknots[level] * dispatch[r,level] for r in 1:N for level in length(turbine(r).powerknots))

        # ------------------------------------------------------------------
        # Flow out
        turbineflow[r=1:N], outflow[r] == sum(turbine(r).flowknots[level] * dispatch[r, level] for level in 1:length(turbine(r).flowknots))

        # Dispatch combination of levels
        dispatched[r=1:N], sum(dispatch[r, level] for level in 1:length(turbine(r).flowknots)) <= 1

        # inflows
        inflows[r=1:N], inflow[r] <= valley_chain[r].inflows[1]
    end)


    # ------------------------------------------------------------------
    #   Objective Function
    stageobjective!(sp, prices[stage, markov_state]*generation_quantity - sum(valley_chain[i].spill_cost * spill[i] for i in 1:N))

end

SDDP.solve(m,
    max_iterations = 10
)
@test isapprox(m.log[end].bound, 846.8, atol=1e-2)
# objective = 846.8
