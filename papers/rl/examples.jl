#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, Gurobi, Ipopt

function example_putterman(; M = 5, N = 10)
    model = SDDP.LinearPolicyGraph(;
        stages = N,
        sense = :Min,
        lower_bound = 0,
        optimizer = Gurobi.Optimizer,
    ) do sp, node
        @variable(sp, x >= 0, SDDP.State, initial_value = M)
        @variable(sp, u >= 0)
        @constraint(sp, x.out == x.in - u)
        @stageobjective(sp, u^2)
        if node == N
            fix(x.out, 0; force = true)
        end
        return
    end
    SDDP.train(model)
    @assert ≈(SDDP.calculate_bound(model), M^2 / N; atol = 1e-6)
    simulations = SDDP.simulate(
        model,
        1,
        [:u],
        sampling_scheme = SDDP.InSampleMonteCarlo(
            max_depth = 10,
            terminate_on_dummy_leaf = false,
        ),
    )
    plt = SDDP.SpaghettiPlot(simulations)
    SDDP.add_spaghetti(plt, title = "u", ymin = 0, ymax = M) do data
        return data[:u]
    end
    SDDP.save(plt)
    return
end

function example_putterman_cyclic(; M = 5, ρ = 0.8)
    model = SDDP.PolicyGraph(
        SDDP.UnicyclicGraph(ρ);
        sense = :Max,
        upper_bound = M / (1 - ρ),
        optimizer = Ipopt.Optimizer,
    ) do sp, node
        @variable(sp, 0 <= x <= M, SDDP.State, initial_value = 0)
        @variable(sp, u >= 0)
        @constraint(sp, x.out == x.in + u)
        @stageobjective(sp, sqrt(1 + u))
        return
    end
    SDDP.train(model)
    simulations = SDDP.simulate(
        model,
        1,
        [:u],
        sampling_scheme = SDDP.InSampleMonteCarlo(
            max_depth = 10,
            terminate_on_dummy_leaf = false,
        ),
    )
    plt = SDDP.SpaghettiPlot(simulations)
    SDDP.add_spaghetti(plt, title = "u", ymin = 0, ymax = M) do data
        return data[:u]
    end
    SDDP.save(plt)
    return
end

function example_tiger_problem()
    ρ = 0.98
    graph = SDDP.Graph(
        :R,
        [:l, :r],
        [(:R => :l, 0.5), (:R => :r, 0.5), (:l => :l, ρ), (:r => :r, ρ)],
    )
    SDDP.add_ambiguity_set(graph, [:l, :r], 1e3)
    model = SDDP.PolicyGraph(
        graph;
        sense = :Min,
        lower_bound = -1000.0,
        optimizer = Gurobi.Optimizer,
    ) do sp, node
        # s: stay, l: open left, r: open right
        @variable(sp, x_s, Bin, SDDP.State, initial_value = 1)
        @variable(sp, x_l, Bin, SDDP.State, initial_value = 0)
        @variable(sp, x_r, Bin, SDDP.State, initial_value = 0)
        @constraint(sp, x_s.out + x_l.out + x_r.out <= 1 - x_l.in - x_r.in)
        @constraint(sp, x_s.out + x_l.out + x_r.out <= x_s.in)
        if node == :l
            @stageobjective(sp, 100 * x_l.in - 10 * x_r.in + x_s.in)
            SDDP.parameterize(sp, [:left, :right], [0.85, 0.15]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        elseif node == :r
            @stageobjective(sp, -10 * x_l.in + 100 * x_r.in + x_s.in)
            SDDP.parameterize(sp, [:left, :right], [0.15, 0.85]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        end
    end
    SDDP.train(
        model;
        iteration_limit = 100,
        log_every_iteration = true,
        cut_deletion_minimum = 10_000,
    )
    simulations = SDDP.simulate(
        model,
        100,
        [:x_l, :x_r, :x_s],
        sampling_scheme = SDDP.InSampleMonteCarlo(
            max_depth = 30,
            terminate_on_dummy_leaf = false,
        ),
    )
    plt = SDDP.SpaghettiPlot(simulations)
    SDDP.add_spaghetti(plt, cumulative = true) do data
        return data[:stage_objective]
    end
    SDDP.add_spaghetti(plt, title = "Stopping state", ymin = 0, ymax = 1) do data
        return data[:x_s].out
    end
    SDDP.add_spaghetti(plt, title = "Open left", ymin = 0, ymax = 1) do data
        return data[:x_l].out
    end
    SDDP.add_spaghetti(plt, title = "Open right", ymin = 0, ymax = 1) do data
        return data[:x_r].out
    end
    SDDP.add_spaghetti(plt, title = "Belief-L", ymin = 0, ymax = 1) do data
        return data[:belief][:l]
    end
    SDDP.save(plt)
end

example_putterman()
example_putterman_cyclic()
example_tiger_problem()
