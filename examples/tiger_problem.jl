#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This example, known in the POMDP literature as the "Tiger Problem," is taken
# from the paper:
#   Kaelbling, L. P., Littmany, M. L., & Cassandra, A. R. (1998). Planning and
#   Acting in Partially Observable Stochastic Domains. Artificial Intelligence,
#   101, 99–134.

# Note: Because of numerical issues associated with the discount factor being
# close to 1, we use Gurobi instead of GLPK. GLPK complains about the basis
# being close to singular.

using SDDP
using Gurobi
using Test

function tiger_problem()
    graph = SDDP.Graph(
        :root_node,
        [:Dl, :Dr, :Hl, :Hr],
        [
            (:root_node => :Dl, 0.5),
            (:root_node => :Dr, 0.5),
            (:Dl => :Hl, 1.0),
            (:Hl => :Dl, 0.98),
            (:Dr => :Hr, 1.0),
            (:Hr => :Dr, 0.98),
        ],
    )
    SDDP.add_ambiguity_set(graph, [:Dl, :Dr], 1e3)
    SDDP.add_ambiguity_set(graph, [:Hl, :Hr], 1e3)

    model = SDDP.PolicyGraph(
        graph,
        lower_bound = -1000.0,
        optimizer = () -> Gurobi.Optimizer(),
    ) do subproblem, node
        set_optimizer_attribute(subproblem, "OutputFlag", 0)
        @variable(subproblem, 0 <= x[[:s, :l, :r]] <= 1, SDDP.State, initial_value = 1)
        if node == :Dl || node == :Dr
            @stageobjective(subproblem, x[:s].out - x[:l].out - x[:r].out)
            @constraints(subproblem, begin
                x[:s].out <= x[:s].in
                x[:l].out + x[:r].out <= x[:s].in
            end)
        elseif node == :Hl
            @stageobjective(subproblem, 100 * x[:l].in - 10 * x[:r].in)
            @constraint(subproblem, x[:s].out <= 1 - x[:l].in - x[:r].in)
            @constraint(subproblem, x[:s].out <= x[:s].in)
            SDDP.parameterize(subproblem, [:left, :right], [0.85, 0.15]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        elseif node == :Hr
            @stageobjective(subproblem, -10 * x[:l].in + 100 * x[:r].in)
            @constraint(subproblem, x[:s].out <= 1 - x[:l].in - x[:r].in)
            @constraint(subproblem, x[:s].out <= x[:s].in)
            SDDP.parameterize(subproblem, [:left, :right], [0.15, 0.85]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        end
        # Dummy constraints to force the state variables to be binary.
        @variable(subproblem, 0 <= u[[:s, :l, :r]] <= 1, Bin)
        @constraint(subproblem, [i = [:s, :l, :r]], x[i].out == u[i])
    end

    # Train the policy.
    SDDP.train(model; iteration_limit = 50, print_level = 1, dashboard = true)

    simulations = SDDP.simulate(
        model,
        100,
        [:x],
        sampling_scheme = SDDP.InSampleMonteCarlo(
            max_depth = 30,
            terminate_on_dummy_leaf = false,
        ),
    )

    plt = SDDP.SpaghettiPlot(simulations)
    SDDP.add_spaghetti(plt, cumulative = true) do data
        data[:stage_objective]
    end
    SDDP.add_spaghetti(plt, title = "Stopping state", ymin = 0, ymax = 1) do data
        data[:x][:s].out
    end
    SDDP.add_spaghetti(plt, title = "Open left", ymin = 0, ymax = 1) do data
        data[:x][:l].out
    end
    SDDP.add_spaghetti(plt, title = "Open right", ymin = 0, ymax = 1) do data
        data[:x][:r].out
    end
    SDDP.add_spaghetti(plt, title = "Belief-L", ymin = 0, ymax = 1) do data
        data[:belief][:Dl] + data[:belief][:Hl]
    end
    SDDP.save(plt)
end

tiger_problem()
