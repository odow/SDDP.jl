#  Copyright 2017-20, Oscar Dowson
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

#=
    The code in this file runs the examples from the paper

    Downward, A., Dowson, O., and Baucke, R. (2018). On the convergence of a
    cutting plane method for multistage stochastic programming problems with
    stagewise dependent price uncertainty. Optimization Online.
=#

using SDDP, GLPK, Test, Gurobi, Plots, StatsPlots

function river_chain_example(;
    ar1::Bool = true,
    N::Int = 2,
    lipschitz = 1e6,
    lower_bound = -50_000,
)
    env = Gurobi.Env()
    model = SDDP.LinearPolicyGraph(
        stages = 12,
        optimizer = () -> Gurobi.Optimizer(env),
        lower_bound = lower_bound,
    ) do sp, t
        set_optimizer_attribute(sp, "OutputFlag", 0)
        flow_knots = [50.0, 60.0, 70.0]
        power_knots = [55.0, 65.0, 70.0]
        b = [
            61.261,
            56.716,
            59.159,
            66.080,
            72.131,
            76.708,
            76.665,
            76.071,
            76.832,
            69.970,
            69.132,
            67.176,
        ]
        Ω = [-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
        @variable(sp, 0 <= volume[1:N] <= 200, SDDP.State, initial_value = 100)
        @variables(sp, begin
            0 <= flow[1:N] <= 70.0
            0 <= spill[1:N]
            0 <= generation
            0 <= dispatch[1:N, 1:3] <= 1
        end)
        @constraints(
            sp,
            begin
                volume[1].out == volume[1].in - flow[1] - spill[1]
                [i = 2:N],
                volume[i].out == volume[i].in - flow[i] - spill[i] + flow[i-1] + spill[i-1]
                generation == sum(power_knots[j] * dispatch[i, j] for i = 1:N, j = 1:3)
                [i = 1:N], flow[i] == sum(flow_knots[j] * dispatch[i, j] for j = 1:3)
                [i = 1:N], sum(dispatch[i, j] for j = 1:3) <= 1
                [i = 1:N], flow[i] <= volume[i].in
            end
        )
        if ar1
            SDDP.add_objective_state(
                sp,
                initial_value = 61.261,
                lower_bound = 40.0,
                upper_bound = 100.0,
                lipschitz = lipschitz,
            ) do p, ω
                if t == 1
                    return p
                else
                    return 0.5 * p[1] + 0.5 * b[t] + ω
                end
            end
            SDDP.parameterize(sp, Ω) do ω
                p′ = SDDP.objective_state(sp)
                @stageobjective(sp, 1_000 * sum(spill) - p′ * generation)
            end
        else
            SDDP.add_objective_state(
                sp,
                initial_value = (61.261, 61.261),
                lower_bound = (40.0, 40.0),
                upper_bound = (100.0, 100.0),
                lipschitz = (lipschitz, lipschitz),
            ) do p, ω
                if t == 1
                    return p
                else
                    return 0.5 * p[1] + 0.5 * b[t] - 0.5 * (p[1] - p[2]) + ω, p[1]
                end
            end
            SDDP.parameterize(sp, Ω) do ω
                p′, p = SDDP.objective_state(sp)
                @stageobjective(sp, 1_000 * sum(spill) - p′ * generation)
            end
        end
    end
    return model
end

function example_one()
    model = river_chain_example(ar1 = true, N = 2)
    SDDP.train(
        model;
        iteration_limit = 2000,
        stopping_rules = [SDDP.Statistical(num_replications = 250, iteration_period = 250)],
    )
    # Now plot the saddle function.
    node = model[6]
    x = 0.0:10.0:200.0
    p = 50.0:5.0:100.0
    A = zeros(length(x), length(p))
    for (j, pj) in enumerate(p)
        for (i, xi) in enumerate(x)
            objective_state_vector =
                SDDP.update_objective_state(node.objective_state, [pj], 0.0)
            subproblem_results = SDDP.solve_subproblem(
                model,
                node,
                Dict(Symbol("volume[1]") => 50.0, Symbol("volume[2]") => xi),
                0.0,
                Tuple{Int,Any}[],
                require_duals = false,
            )
            A[i, j] = subproblem_results.objective
        end
    end
    open("saddlefunction.dat", "w") do io
        for (i, xi) in enumerate(x)
            for (j, pj) in enumerate(p)
                println(io, xi, " ", pj, " ", A[i, j])
            end
            println(io)
        end
    end
end

function example_two()
    model = river_chain_example(ar1 = false, N = 5, lower_bound = -150_000.0)
    SDDP.train(model; iteration_limit = 5_000)
    simulations = SDDP.simulate(model, 1_000, [:generation, :volume])

    profits = map(sim -> sum(s[:stage_objective] for s in sim), simulations)

    min_profit, min_index = findmin(profits)
    max_profit, max_index = findmax(profits)

    function pretty_plot(f, simulations; kwargs...)
        plt = SDDP.publication_plot(f, simulations)
        high = f.(simulations[max_index])
        plot!(high, linestyle = :dash, color = :red, width = 3)
        low = f.(simulations[min_index])
        plot!(low, linestyle = :dot, color = :green, width = 3)
        plot!(; legend = false, xlabel = "Stage", kwargs...)
        return plt
    end

    spot_price = pretty_plot(
        simulations;
        title = "(a) Spot Price",
        ylabel = "Spot Price\n(\$/MWh)",
    ) do data
        return data[:objective_state][1]
    end

    stored_energy = pretty_plot(
        simulations;
        title = "(b) Total Stored Energy",
        ylabel = "Volume\n(m³)",
    ) do data
        return sum(data[:volume][i].out for i = 1:5)
    end

    total_generation = pretty_plot(
        simulations;
        title = "(b) Total Generation",
        ylabel = "Energy\n(MWh)",
    ) do data
        return data[:generation]
    end

    profit = StatsPlots.density(
        profits ./ 1_000;
        title = "(d) Distribution of Profit",
        xlabel = "Profit (000's)",
        width = 3,
        color = "#00467F",
        legend = false,
        yticks = false,
    )
    Plots.vline!([min_profit / 1_000], color = :green, linestyle = :dot, width = 3)
    Plots.vline!([max_profit / 1_000], color = :red, linestyle = :dash, width = 3)
    plot(spot_price, stored_energy, total_generation, profit)
    savefig("example_two.pdf")
end

function example_lipschitz()
    A = zeros(500, 5)
    for (col, α) in enumerate([0.0, 10.0, 100.0, 1_000])
        model = river_chain_example(ar1 = false, N = 5, lipschitz = α)
        SDDP.train(model; iteration_limit = 500)
        log = model.most_recent_training_results.log
        for row = 1:500
            A[row, col] = log[row].bound
        end
    end
    open("lipschitz_experiment.dat", "w") do io
        for row = 1:500
            print(rpad(row, 4))
            for col = 1:5
                print(A[row, col], " ")
            end
            println()
        end
    end
end

if length(ARGS) > 0
    if ARGS[1] == "example_one"
        example_one()
    elseif ARGS[1] == "example_two"
        example_two()
    elseif ARGS[1] == "example_lipschitz"
        example_lipschitz()
    end
end
