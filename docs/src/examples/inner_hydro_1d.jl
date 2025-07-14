#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Hydro-thermal with inner approximation

# This is a simple version of the hydro-thermal scheduling problem. The goal is
# to operate one hydro-dam and two thermal plants over time in the face of
# inflow uncertainty.

# The purpose of this tutorial is to provide a demonstration of the experimental
# `SDDP.Inner` submodule.
#
# !!! warning
#     The `SDDP.Inner` code in this example is experimental and the API may
#     change in any future release.

using SDDP
using Test

import HiGHS
import Random

function test_inner_hydro_1d()
    Random.seed!(2)
    stages = 4
    Ω = max.(40 .+ 20.0 * randn(10), 0.0)
    risk_measure = SDDP.EAVaR(; lambda = 0.5, beta = 0.2)
    function build_subproblem(sp::JuMP.Model, t::Int)
        @variable(sp, 0 <= x_vol <= 100, SDDP.State, initial_value = 83.222)
        @variable(sp, 0 <= u_gh <= 60)
        @variable(sp, 0 <= u_spill <= 200)
        @variable(sp, 0 <= u_gt[1:2] <= 15)
        @variable(sp, 0 <= u_def <= 75)
        @variable(sp, w_inflow == 0.0)
        @constraint(sp, sum(u_gt) + u_gh + u_def == 75)
        @constraint(sp, x_vol.in + w_inflow - u_gh - u_spill == x_vol.out)
        if t > 1
            SDDP.parameterize(w -> JuMP.fix(w_inflow, w), sp, Ω)
        end
        @stageobjective(sp, u_spill + 5 * u_gt[1] + 10 * u_gt[2] + 50 * u_def)
    end
    println("Building and solving primal outer model for lower bounds")
    model = SDDP.LinearPolicyGraph(
        build_subproblem;
        stages,
        sense = :Min,
        optimizer = HiGHS.Optimizer,
        lower_bound = 0.0,
    )
    SDDP.train(model; iteration_limit = 10, risk_measure)
    lower_bound = SDDP.calculate_bound(model)

    simulations =
        SDDP.simulate(model, 500, [:x_vol, :u_gh, :u_spill, :u_gt, :u_def])
    objs = [sum(data[:stage_objective] for data in sim) for sim in simulations]
    μ, ci = round.(SDDP.confidence_interval(objs); digits = 2)
    println("Building and solving inner model for upper bounds:")
    inner_model, upper_bound = SDDP.Inner.inner_dp(
        build_subproblem,
        model;
        stages,
        sense = :Min,
        optimizer = HiGHS.Optimizer,
        lower_bound = 0.0,
        risk_measure,
        bellman_function = SDDP.Inner.InnerBellmanFunction(
            t -> 50.0 * (stages - t);
            upper_bound = t -> 75.0 * 50.0 * (stages - t),
            vertex_type = SDDP.SINGLE_CUT,
        ),
    )
    @test lower_bound <= upper_bound
    println()
    println("Bounds:")
    println("  Risk-neutral confidence interval: ", μ, " ± ", ci)
    println("  Risk-adjusted lower bound: ", round(lower_bound; digits = 2))
    println("  Risk-adjusted upper bound: ", round(upper_bound; digits = 2))
    return
end

test_inner_hydro_1d()
