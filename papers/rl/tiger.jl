#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP
import Gurobi
import Plots
import Random

function example_tiger_problem(ε::Float64 = 0.35; create_plot::Bool = true)
    ρ = 0.98
    graph = SDDP.Graph(
        :R,
        [:l, :r],
        [(:R => :l, 0.5), (:R => :r, 0.5), (:l => :l, ρ), (:r => :r, ρ)],
    )
    SDDP.add_ambiguity_set(graph, [:l, :r], 1e4)
    model = SDDP.PolicyGraph(
        graph;
        sense = :Min,
        lower_bound = -10.0,
        optimizer = Gurobi.Optimizer,
    ) do sp, node
        set_attribute(sp, "IntegralityFocus", 1)
        # s: stay, l: open left, r: open right
        @variable(sp, x_s, Bin, SDDP.State, initial_value = 1)
        @variable(sp, x_l, Bin, SDDP.State, initial_value = 0)
        @variable(sp, x_r, Bin, SDDP.State, initial_value = 0)
        @constraint(sp, x_s.out + x_l.out + x_r.out <= 1 - x_l.in - x_r.in)
        @constraint(sp, x_s.out + x_l.out + x_r.out <= x_s.in)
        if node == :l
            @stageobjective(sp, 100 * x_l.in - 10 * x_r.in + x_s.in)
            SDDP.parameterize(sp, [:left, :right], [0.5 + ε, 0.5 - ε]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        elseif node == :r
            @stageobjective(sp, -10 * x_l.in + 100 * x_r.in + x_s.in)
            SDDP.parameterize(sp, [:left, :right], [0.5 - ε, 0.5 + ε]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        end
    end
    Random.seed!(123)
    SDDP.train(
        model;
        iteration_limit = 100,
        log_every_iteration = true,
        # cut_deletion_minimum = 1_000,
        duality_handler = SDDP.LagrangianDuality(),
        # duality_handler = SDDP.BanditDuality(
        #     SDDP.ContinuousConicDuality(),
        #     SDDP.LagrangianDuality(),
        # ),
    )
    lower_bound = SDDP.calculate_bound(model)
    Random.seed!(456)
    sampling_scheme = SDDP.InSampleMonteCarlo(
        max_depth = 30,
        terminate_on_dummy_leaf = false,
    )
    simulations = SDDP.simulate(model, 100, [:x_l, :x_r]; sampling_scheme);
    objectives = map(simulations) do simulation
        return sum(d[:stage_objective] for d in simulation)
    end
    μ, σ = SDDP.confidence_interval(objectives)
    if !create_plot
        println("lower_bound = $lower_bound")
        println("upper_bound = $μ ± $σ")
        return lower_bound, μ, σ
    end
    simulations = simulations[1:100]
    belief_plot = Plots.plot(;
        xlabel = "Time step",
        ylabel = "Belief(Left)",
        legend = false,
        ymajorgrid = true,
    )
    plot = Plots.plot(;
        xlabel = "Time step",
        ylabel = "# hear left - hear right",
        legend = false,
        ymajorgrid = true,
    )
    for simulation in simulations
        b = Float64[0.5]
        y = Int[0]
        for d in simulation
            push!(b, d[:belief][:l])
            if d[:noise_term] == :left
                push!(y, y[end] + 1)
            else
                push!(y, y[end] - 1)
            end
            if d[:x_l].out > 0.5 || d[:x_r].out > 0.5
                break
            end
        end
        Plots.plot!(
            belief_plot,
            0:length(b)-1,
            b;
            color = :grey,
            linewidth = 3,
            alpha = 0.2,
        )
        Plots.plot!(
            plot,
            0:length(y)-1,
            y;
            color = :grey,
            linewidth = 3,
            alpha = 0.2,
        )
        if (i = findfirst(d -> d[:x_l].out > 0.5, simulation)) !== nothing
            Plots.scatter!([i], [y[i+1]], color = "#43a047", markersize = 6)
        end
        if (i = findfirst(d -> d[:x_r].out > 0.5, simulation)) !== nothing
            Plots.scatter!(
                [i],
                [y[i+1]];
                marker = :x,
                markersize = 6,
                markerstrokewidth = 3,
                color = "#e53935",
            )
        end
    end
    Plots.plot(belief_plot, plot, layout = (2, 1), dpi = 400)
    Plots.savefig("tiger_problem_$ε.pdf")
    return lower_bound, μ, σ
end

# data = [
#     example_tiger_problem(ε; create_plot = false) for ε in [0.35, 0.4, 0.45]
# ]

d_35 = example_tiger_problem(0.35; create_plot = true)
