# Copyright (c) 2017-26: Oscar Dowson and SDDP.jl contributors.          #src
#                                                                        #src
# This Source Code Form is subject to the terms of the Mozilla Public    #src
# License, v2.0. If a copy of the MPL was not distributed with this file #src
# You can obtain one at http://mozilla.org/MPL/2.0/.                     #src

# # Example: production planning

using SDDP
import HiGHS
import Plots

P, M = ["Seattle", "San-Diego"], ["New-York", "Chicago", "Topeka"]
Ω = [
    Dict("New-York" => 300, "Chicago" => 300, "Topeka" => 300),
    Dict("New-York" => 350, "Chicago" => 320, "Topeka" => 310),
    Dict("New-York" => 320, "Chicago" => 390, "Topeka" => 350),
    Dict("New-York" => 250, "Chicago" => 220, "Topeka" => 330),
    Dict("New-York" => 290, "Chicago" => 200, "Topeka" => 290),
]
c_capacity = Dict("Seattle" => 350, "San-Diego" => 600)
c_ship = Containers.DenseAxisArray([2.6 1.7 1.8; 2.5 1.8 1.4], P, M)
model = SDDP.LinearPolicyGraph(;
    stages = 10,
    sense = :Min,
    optimizer = HiGHS.Optimizer,
    lower_bound = 0.0,
) do sp, t
    ## State variable: inventory levels
    @variable(sp, x[p in P] >= 0, SDDP.State, initial_value = 0)
    ## Control variable: quantity to produce at each plant
    @variable(sp, 0 <= u_prod[p in P] <= c_capacity[p])
    ## Control variable: quantity to ship from plant p to market m
    @variable(sp, u_ship[p in P, m in M] >= 0)
    ## Control variable: unmet demand
    @variable(sp, u_unmet[m in M] >= 0)
    ## Random variable: demand in each market
    @variable(sp, w[m in M] == 0)
    if t > 1
        ## In the first-stage there is no uncertainty
        SDDP.parameterize(sp, Ω) do ω
            for m in M
                fix(w[m], ω[m])
            end
            return
        end
    end
    ## Constraint: can ship only from existing inventory
    @constraint(sp, [p in P], sum(u_ship[p, :]) <= x[p].in)
    ## Constraint: balance inventory at each plant
    @constraint(
        sp,
        [p in P],
        x[p].out == x[p].in + u_prod[p] - sum(u_ship[p, :]),
    )
    ## Constraint: balance demand in each market
    @constraint(sp, [m in M], sum(u_ship[:, m]) == w[m] - u_unmet[m])
    @stageobjective(
        sp,
        ## Cost of production and holding inventory
        sum(1.0 * u_prod[p] + 0.01 * x[p].out for p in P) +
        ## Cost of unmet demand
        sum(10 * u_unmet[m] for m in M) +
        ## Cost of shipping
        sum(c_ship[p, m] * u_ship[p, m] for p in P, m in M),
    )
end

#-

SDDP.train(model; risk_measure = SDDP.Expectation())

#-

simulations = SDDP.simulate(model, 50, [:x, :u_prod]);

#-

Plots.plot(
    SDDP.publication_plot(
        simulations;
        title = "San-Diego",
        ylabel = "u_prod",
    ) do data
        return data[:u_prod]["San-Diego"]
    end,
    SDDP.publication_plot(simulations; title = "Seattle") do data
        return data[:u_prod]["Seattle"]
    end,
    SDDP.publication_plot(simulations; ylabel = "x", xlabel = "Week") do data
        return data[:x]["San-Diego"].out
    end,
    SDDP.publication_plot(simulations; xlabel = "Week") do data
        return data[:x]["Seattle"].out
    end,
)
