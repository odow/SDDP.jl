#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP
import Gurobi
import HiGHS
import Plots
import Random
import StatsPlots

mutable struct StoppingForwardPass <: SDDP.AbstractForwardPass
    name::Symbol
end

function SDDP.forward_pass(
    model::SDDP.PolicyGraph,
    options::SDDP.Options,
    f_pass::StoppingForwardPass,
)
    pass = SDDP.forward_pass(model, options, SDDP.DefaultForwardPass())
    index = findfirst(s -> s[f_pass.name] < 0.5, pass.sampled_states)
    n = 3 + something(index, length(pass.sampled_states))
    subset(x, n) = x[1:min(length(x), n)]
    return (
        scenario_path = subset(pass.scenario_path, n),
        sampled_states = subset(pass.sampled_states, n),
        objective_states = pass.objective_states,
        belief_states = subset(pass.belief_states, n),
        cumulative_value = pass.cumulative_value,
    )
end

function example_tiger_problem(ε::Float64; create_plot::Bool)
    ρ = 0.95
    graph = SDDP.Graph(
        :R,
        [:l, :r],
        [(:R => :l, 0.5), (:R => :r, 0.5), (:l => :l, ρ), (:r => :r, ρ)],
    )
    if ε < 0.5
        SDDP.add_ambiguity_set(graph, [:l, :r], 1e3)
    end
    model = SDDP.PolicyGraph(
        graph;
        sense = :Min,
        lower_bound = -10.0,
        optimizer = Gurobi.Optimizer,
    ) do sp, node
        # s: stay, l: open left, r: open right
        @variable(sp, x_s, Bin, SDDP.State, initial_value = 1)
        @variable(sp, x_l, Bin, SDDP.State, initial_value = 0)
        @variable(sp, x_r, Bin, SDDP.State, initial_value = 0)
        @constraint(sp, x_s.out + x_l.out + x_r.out <= 1 - x_l.in - x_r.in)
        @constraint(sp, x_s.out + x_l.out + x_r.out == x_s.in)
        if node == :l
            @stageobjective(sp, 100 * x_l.in - 10 * x_r.in + x_s.in)
            SDDP.parameterize(sp, [:left, :right], [0.5+ε, 0.5-ε]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        elseif node == :r
            @stageobjective(sp, -10 * x_l.in + 100 * x_r.in + x_s.in)
            SDDP.parameterize(sp, [:left, :right], [0.5-ε, 0.5+ε]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        end
    end
    Random.seed!(1234)
    SDDP.train(
        model;
        iteration_limit = ε < 0.15 ? 200 : 100,
        log_every_iteration = true,
        cut_deletion_minimum = 1_000,
        duality_handler = SDDP.LagrangianDuality(),
        forward_pass = StoppingForwardPass(:x_s),
    )
    lower_bound = SDDP.calculate_bound(model)
    Random.seed!(456)
    sampling_scheme = SDDP.InSampleMonteCarlo(
        max_depth = 50,
        terminate_on_dummy_leaf = false,
    )
    simulations =
        SDDP.simulate(model, 1_000, [:x_s, :x_l, :x_r]; sampling_scheme)
    objectives = map(simulations) do simulation
        return sum(
            ρ^(t-1) * d[:stage_objective] for (t, d) in enumerate(simulation)
        )
    end
    if !create_plot
        μ, σ = SDDP.confidence_interval(objectives)
        println("lower_bound = $lower_bound")
        println("upper_bound = $μ ± $σ")
        return lower_bound, objectives, model
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
    return lower_bound, objectives, model
end

function run_experiment()
    data = Dict()
    for ε in [0.1, 0.2, 0.3, 0.4, 0.5]
        data[ε] = example_tiger_problem(ε; create_plot = false)
    end

    x = sort(collect(keys(data)))
    lb = [data[xi][1] for xi in x]
    μ = [data[xi][2] for xi in x]
    box_y = reduce(vcat, μ)
    box_x = reduce(vcat, [fill(0.5 - x[i], length(μ[i])) for i in 1:length(x)])

    StatsPlots.violin(
        box_x[box_y .<= 60],
        box_y[box_y .<= 60] .+ 0.1 * rand(sum(box_y .<= 60));
        bar_width = 0.05,
        ylims = (-10, 25),
        xlabel = "False positive rate",
        ylabel = "Objective value",
        label = false,
        color = :grey,
        alpha = 0.5,
    )
    Plots.scatter!(
        0.5 .- x,
        lb;
        label = "Lower bound",
        color = :black,
        linewidth = 3,
    )
    Plots.savefig("tiger_problem_violin.pdf")
    return
end

run_experiment()
