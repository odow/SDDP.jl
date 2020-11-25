#  Copyright 2019-21, Oscar Dowson, Lingquan Ding (@lingquant).
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This example is based on one from MSPPy:
# https://github.com/lingquant/msppy/blob/dc85a2e8fa5243b3d5096d59085d9caad3ff2ede/examples/hydro_thermal/julia/test.jl
#
# The original author was Lingquan Ding (@lingquant), but it was modified by
# Oscar Dowson (@odow) to meet the latest SDDP.jl syntax.
#
# The original model and data is from:
#
# Shapiro, A., Tekaya, W., da Costa, J. P., & Soares, M. P. (2013). Risk neutral
# and risk averse stochastic dual dynamic programming method. European journal
# of operational research, 224(2), 375–391.

using SDDP

import DelimitedFiles
import Gurobi
import JSON
import Plots
import Random
import Statistics
import StatsPlots

function _build_model(data::Dict)
    N_THERMAL = length.(data["thermal_obj"])
    @assert N_THERMAL == length.(data["thermal_lb"])
    @assert N_THERMAL == length.(data["thermal_ub"])
    I = [1, 2, 3, 4]
    K(i) = 1:N_THERMAL[i]
    IM = 5
    I_IM = union(I, IM)
    env = Gurobi.Env()
    function gurobi_optimizer()
        model = Gurobi.Optimizer(env)
        MOI.set(model, MOI.Silent(), true)
        return model
    end
    model = SDDP.LinearPolicyGraph(
        stages = 60,
        lower_bound = 0.0,
        optimizer = gurobi_optimizer,
    ) do sp, t
        set_silent(sp)
        month = t % 12 == 0 ? 12 : t % 12
        @variable(
            sp,
            0 <= v[i in I] <= data["storedEnergy_ub"][i] / 1_000,
            SDDP.State,
            initial_value = data["storedEnergy_initial"][i] / 1_000,
        )
        @variable(sp, s[i in I] >= 0.0)
        @variable(sp, 0.0 <= q[i in I] <= data["hydro_ub"][i] / 1_000)
        @variable(
            sp,
            g[i in I, k in K(i)],
            lower_bound = data["thermal_lb"][i][k] / 1_000,
            upper_bound = data["thermal_ub"][i][k] / 1_000,
        )
        @variable(
            sp,
            0 <= ex[i in I_IM, j in I_IM] <= data["exchange_ub"][i][j] / 1_000,
        )
        @variable(
            sp,
            df[i in I, j in I] >= 0,
            upper_bound =
                data["demand"][month][i] * data["deficit_ub"][j] / 1_000,
        )
        @stageobjective(
            sp,
            sum(
                data["deficit_obj"][i] * sum(df[i, :]) +
                sum(data["thermal_obj"][i][k] * g[i, k] for k in K(i)) for
                i in I
            )
        )
        @constraint(
            sp,
            [i in I],
            q[i] + sum(g[i, k] for k in K(i)) + sum(df[i, :]) + sum(ex[:, i]) - sum(ex[i, :]) ==
            data["demand"][month][i] / 1_000
        )
        @constraint(
            sp,
            balance[i in I],
            v[i].out ==
            v[i].in + data["inflow_initial"][i] / 1_000 - s[i] - q[i]
        )
        @constraint(sp, sum(ex[:, IM]) == sum(ex[IM, :]))
        if t == 1
            # Deterministic first stage with `inflow_initial`.
        else
            r = (t - 1) % 12 == 0 ? 12 : (t - 1) % 12
            Ω = data["scenarios"]
            # To simplify things. Don't use all the scenarios!
            num_scenarios = min(40, length(Ω[1][r]))
            SDDP.parameterize(sp, 1:num_scenarios) do ω
                for i in 1:4
                    set_normalized_rhs(balance[i], Ω[i][r][ω] / 1_000)
                end
            end
        end
    end
    return model
end

function _train_model(model::SDDP.PolicyGraph, gamma::Float64)
    risk_measure = isfinite(gamma) ? SDDP.Entropic(gamma) : SDDP.WorstCase()
    # Set the same random seed for all runs!
    Random.seed!(1234)
    psr_sampling_scheme = SDDP.PSRSamplingScheme(1_000)
    try
        SDDP.train(
            model;
            iteration_limit = 10_000,
            risk_measure = risk_measure,
            sampling_scheme = psr_sampling_scheme,
            log_file = "$(gamma).log",
        )
    catch ex
        @info "Ignorring error"
        println("$(ex)")
    end
    # Set the same random seed for all runs!
    Random.seed!(12345)
    simulations = SDDP.simulate(
        model,
        1_000,
        [:v, :df, :g];
        sampling_scheme = psr_sampling_scheme,
    )
    stage_objectives = [
        round(simulations[i][t][:stage_objective]; digits = 1) for
        i in 1:length(simulations) for t in 1:length(simulations[i])
    ]
    stage_objectives =
        reshape(stage_objectives, length(simulations[1]), length(simulations))
    open("stage_objectives_$(gamma).csv", "w") do io
        return DelimitedFiles.writedlm(
            io,
            collect(transpose(stage_objectives)),
            ',',
        )
    end
    return
end

function _plot_cumulative_density(
    filenames = [
        "stage_objectives_0.000000.csv",
        "stage_objectives_0.000010.csv",
        "stage_objectives_5.0e-5.csv",
        "stage_objectives_0.000100.csv",
    ],
)
    Plots.plot()
    line_styles = [:solid, :dash, :dot, :dashdot]
    for (i, f) in enumerate(filenames)
        X = sum(DelimitedFiles.readdlm(f, ','); dims = 2)[:]
        x = parse(Float64, match(r"stage\_objectives\_(.+)\.csv", f)[1])
        x = if x ≈ 0.0
            "0"
        elseif x ≈ 5e-5
            "5 \\times 10^{-5}"
        else
            "1 \\times 10^{$(Int(log10(x)))}"
        end
        StatsPlots.cdensity!(
            X,
            label = "\$\\gamma = $(x)\$",
            style = line_styles[i],
            color = :black,
        )
    end
    p = Plots.plot!(
        xlabel = "Cost [\$]",
        ylabel = "Cumulative density",
        legend = :bottomright,
        size = (400, 300),
    )
    Plots.savefig("cumulative_density.pdf")
    return p
end

function _plot_objectives(filename::String)
    quantiles = [0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0]
    matrix = DelimitedFiles.readdlm(filename, ',')
    A = mapreduce(
        i -> Statistics.quantile(matrix[:, i], quantiles)',
        vcat,
        1:size(matrix, 2),
    )
    x = parse(Float64, match(r"stage\_objectives\_(.+)\.csv", filename)[1])
    x = if x ≈ 0.0
        "0"
    elseif x ≈ 5e-5
        "5 \\times 10^{-5}"
    else
        "10^{$(Int(log10(x)))}"
    end
    Plots.plot(
        A[:, 4],
        ribbon = (A[:, 4] .- A[:, 1], A[:, 7] .- A[:, 4]),
        color = :black,
        fillalpha = 0.2,
        legend = false,
        title = "\$\\gamma = $(x)\$",
        size = (300, 400),
    )
    Plots.plot!(
        A[:, 4],
        ribbon = (A[:, 4] .- A[:, 2], A[:, 6] .- A[:, 4]),
        color = :black,
        fillalpha = 0.2,
    )
    Plots.plot!(
        A[:, 4],
        ribbon = (A[:, 4] .- A[:, 3], A[:, 5] .- A[:, 4]),
        color = :black,
        fillalpha = 0.2,
    )
    return Plots.plot!(
        A[:, 4],
        color = :black,
        xlabel = "Stages",
        ylabel = "Stage objective (\$)",
        ylim = (0, 4e4),
    )
end

function _plot_objectives(
    filenames::Vector{String} = [
        "stage_objectives_0.000000.csv",
        "stage_objectives_0.000010.csv",
        "stage_objectives_5.0e-5.csv",
        "stage_objectives_0.000100.csv",
    ],
)
    p = Plots.plot(
        _plot_objectives.(filenames)...,
        size = (1200, 800),
        margin = 5Plots.mm,
    )
    Plots.savefig("objectives.pdf")
    return p
end

function _print_help()
    return println(
        """
usage: julia [-p N] [--project=.] model.jl [--gamma=<value>] [--help] [--plot]

Solve the hydro-thermal scheduling problem with the Entropic risk measure
parameterized with γ = <value>.

Use `-p N` to run the SDDP solver in parallel mode over `N` processors.

Use `--project=.` to reproduce the example in the paper using the provided
`Manifest.toml` and `Project.toml` files.

Examples:

    julia model.jl --gamma=0.5
    julia -p 4 --project=. model.jl --gamma=0.5
    julia -p 4 --project=. model.jl --gamma=5.0
    julia --project=. model.jl --plot
""",
    )
end

function main(args)
    if length(args) == 0
        _print_help()
    elseif findfirst(isequal("--help"), args) !== nothing
        _print_help()
    elseif findfirst(isequal("-h"), args) !== nothing
        _print_help()
    else
        i = findfirst(arg -> startswith(arg, "--gamma="), args)
        if i !== nothing
            data = JSON.parsefile(joinpath(@__DIR__, "data.json"))
            model = _build_model(data)
            gamma = parse(Float64, replace(args[i], "--gamma=" => ""))::Float64
            _train_model(model, gamma)
        end
        if findfirst(isequal("--plot"), args) !== nothing
            _plot_cumulative_density()
            _plot_objectives()
        end
    end
end

# ============================================================================ #
#                              Script entry point!                             #
# ============================================================================ #

main(ARGS)
