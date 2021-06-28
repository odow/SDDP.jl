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
#
# This code is used to build the figures contained in the paper
#
# On solving multistage stochastic programs with the entropic risk measure.
# Dowson, O., Morton, D.P., and Pagnoncelli, B.K.

using SDDP

import DelimitedFiles
import Gurobi
import JSON
import Plots
import Random
import Statistics

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
                sum(data["thermal_obj"][i][k] * g[i, k] for k in K(i))
                for i in I
            )
        )
        @constraint(
            sp,
            [i in I],
            q[i] + sum(g[i, k] for k in K(i)) + sum(df[i, :]) + sum(ex[:, i]) -
            sum(ex[i, :]) == data["demand"][month][i] / 1_000
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
            num_scenarios = length(Ω[1][r])
            SDDP.parameterize(sp, 1:10) do ω
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
    try
        SDDP.train(
            model;
            iteration_limit = 1_000, # 5_000,
            risk_measure = risk_measure,
            forward_pass = SDDP.RiskAdjustedForwardPass(
                forward_pass = SDDP.DefaultForwardPass(),
                risk_measure = SDDP.Entropic(gamma),
                resampling_probability = 0.1,
                rejection_count = 5,
            ),
            log_file = "$(gamma).log",
        )
    catch ex
        @info "Ignorring error"
        println("$(ex)")
    end
    # Set the same random seed for all runs!
    Random.seed!(12345)
    simulations = SDDP.simulate(model, 2_000, [:v, :df, :g])
    stage_objectives = [
        round(simulations[i][t][:stage_objective]; digits = 1)
        for i = 1:length(simulations)
        for t = 1:length(simulations[i])
    ]
    stage_objectives = reshape(
        stage_objectives,
        length(simulations[1]),
        length(simulations),
    )
    open("stage_objectives_$(gamma).csv", "w") do io
        DelimitedFiles.writedlm(io, collect(transpose(stage_objectives)), ',')
    end
    return
end

function _plot_results()
    quantiles = [0.01, 0.1, 0.5, 0.9, 0.99]
    filenames = filter(f -> endswith(f, ".csv"), readdir(@__DIR__))
    per_stage = Dict{Float64, Matrix{Float64}}()
    total = Dict{Float64, Vector{Float64}}()
    for file in filenames
        key = parse(Float64, match(r"stage\_objectives\_(.+)\.csv", file)[1])
        matrix = DelimitedFiles.readdlm(file, ',')
        per_stage[key] = mapreduce(
            i -> Statistics.quantile(matrix[:, i], quantiles)',
            vcat,
            1:size(matrix, 2),
        )
        total[key] = Statistics.quantile(sum(matrix; dims = 2)[:], quantiles)
    end
    for (i, q) in enumerate(quantiles)
        Plots.plot()
        for (γ, X) in per_stage
            Plots.plot!(X[:, i], label = "$γ")
        end
        Plots.savefig("quantile_$(q).pdf")
    end
    # x = sort(collect(keys(total)))
    # y = mapreduce(xi -> total[xi]', vcat, x)
    # Plots.plot(x, y, labels = false, xscale=:log10)
    # Plots.savefig("end-of-horizon.pdf")
    return
end

function _print_help()
    println("""
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
    """)
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
            _plot_results()
        end
    end
end

# ============================================================================ #
#                              Script entry point!                             #
# ============================================================================ #

main(ARGS)
