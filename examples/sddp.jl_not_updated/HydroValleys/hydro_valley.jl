#  Copyright 2017,  Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
This problem is a version of the hydro-thermal scheduling problem. The goal is
to operate two hydro-dams in a valley chain over time in the face of inflow
and price uncertainty.

Turbine response curves are modelled by piecewise linear functions which map the
flow rate into a power. These can be controlled by specifying the breakpoints
in the piecewise linear function as the knots in the Turbine struct.

The model can be created using the hydrovalleymodel function. It has a few
keyword arguments to allow automated testing of the library.
`hasstagewiseinflows` determines if the RHS noise constraint should be added.
`hasmarkovprice` determines if the price uncertainty (modelled by a markov
chain) should be added.

In the third stage, the markov chain has some unreachable states to test
some code-paths in the library.

We can also set the sense to :Min or :Max (the objective and bound are
flipped appropriately).
=#
using SDDP, JuMP, Clp

struct Turbine
    flowknots::Vector{Float64}
    powerknots::Vector{Float64}
end

struct Reservoir
    min::Float64
    max::Float64
    initial::Float64
    turbine::Turbine
    spill_cost::Float64
    inflows::Vector{Float64}
end

function hydrovalleymodel(;
        hasstagewiseinflows::Bool=true,
        hasmarkovprice::Bool=true,
        sense::Symbol=:Max
    )

    valley_chain = [
        Reservoir(0, 200, 200, Turbine([50, 60, 70], [55, 65, 70]), 1000, [0, 20, 50]),
        Reservoir(0, 200, 200, Turbine([50, 60, 70], [55, 65, 70]), 1000, [0, 0,  20])
    ]

    turbine(i) = valley_chain[i].turbine

    # Prices[t, markov state]
    prices = [
        1 2 0;
        2 1 0;
        3 4 0
    ]

    # Transition matrix
    if hasmarkovprice
        transition = Array{Float64, 2}[
            [ 1.0 ]',
            [ 0.6 0.4 ],
            [ 0.6 0.4 0.0; 0.3 0.7 0.0]
        ]
    else
        transition = [ones(Float64, (1,1)) for t in 1:3]
    end

    flipobj = (sense == :Max) ? 1.0 : -1.0
    lower = (sense == :Max) ? -Inf : -1e6
    upper = (sense == :Max) ? 1e6 : Inf

    N = length(valley_chain)

    # Initialise SDDP Model
    m = SDDP.MarkovianPolicyGraph(
            sense                 = sense
            , lower_bound         = lower
            , upper_bound         = upper
            , transition_matrices = transition
            , optimizer           = with_optimizer(Clp.Optimizer)) do subproblem, node
        t, markov_state = node

        # ------------------------------------------------------------------
        #   SDDP State Variables
        # Level of upper reservoir
        @variable(subproblem, valley_chain[r].min <= reservoir[r=1:N] <= valley_chain[r].max, SDDP.State, initial_value = valley_chain[r].initial)

        # ------------------------------------------------------------------
        #   Additional variables
        @variables(subproblem, begin
            outflow[r=1:N]      >= 0
            spill[r=1:N]        >= 0
            inflow[r=1:N]       >= 0
            generation_quantity >= 0 # Total quantity of water
            # Proportion of levels to dispatch on
            0 <= dispatch[r=1:N, level=1:length(turbine(r).flowknots)] <= 1
            rainfall[i=1:N]
        end)

        # ------------------------------------------------------------------
        # Constraints
        @constraints(subproblem, begin
            # flow from upper reservoir
            reservoir[1].out == reservoir[1].in + inflow[1] - outflow[1] - spill[1]

            # other flows
            flow[i=2:N], reservoir[i].out == reservoir[i].in + inflow[i] - outflow[i] - spill[i] + outflow[i-1] + spill[i-1]

            # Total quantity generated
            generation_quantity == sum(turbine(r).powerknots[level] * dispatch[r,level] for r in 1:N for level in 1:length(turbine(r).powerknots))

            # ------------------------------------------------------------------
            # Flow out
            turbineflow[r=1:N], outflow[r] == sum(turbine(r).flowknots[level] * dispatch[r, level] for level in 1:length(turbine(r).flowknots))

            # Dispatch combination of levels
            dispatched[r=1:N], sum(dispatch[r, level] for level in 1:length(turbine(r).flowknots)) <= 1
        end)

        # rainfall noises
        if hasstagewiseinflows && t > 1 # in future stages random inflows
            @constraint(subproblem, inflow_noise[i=1:N], inflow[i] <= rainfall[i])

            SDDP.parameterize(subproblem, [(valley_chain[1].inflows[i], valley_chain[2].inflows[i]) for i = 1:length(transition)]) do ω
                for i in 1:N
                    JuMP.fix(rainfall[i], ω[i])
                end
            end
        else # in the first stage deterministic inflow
            @constraint(subproblem, initial_inflow_noise[i=1:N], inflow[i] <= valley_chain[i].inflows[1])
        end

        # ------------------------------------------------------------------
        #   Objective Function
        if hasmarkovprice
            @stageobjective(subproblem, flipobj * (prices[t, markov_state] * generation_quantity - sum(valley_chain[i].spill_cost * spill[i] for i in 1:N)))
        else
            @stageobjective(subproblem, flipobj * (prices[t, 1] * generation_quantity - sum(valley_chain[i].spill_cost * spill[i] for i in 1:N)))
        end
    end
end
