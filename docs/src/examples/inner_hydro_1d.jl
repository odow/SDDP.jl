#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Hydro 1d

# This is a simple version of the hydro-thermal scheduling problem.
# The goal is # to operate one hydro-dam and two thermal plants over time
# in the face of inflow uncertainty.

using SDDP, HiGHS, Test, Random, Statistics
import SDDP: Inner

function test_inner_hydro_1d()
    ## Parameters
    ## Model
    nstages = 4
    inivol = 83.222

    ## Uncertainty
    nscen = 10
    Random.seed!(2)
    inflows = max.(40 .+ 20 * randn(nscen), 0)

    ## Risk aversion
    ## EAVaR(;lambda=1.0, beta=1.0) corresponds to
    ##   ρ = λ * E[x] + (1 - λ) * AV@R(β)[x]
    ## and β = 1.0 makes AV@R = E
    lambda = 0.5
    beta = 0.20
    rho = SDDP.EAVaR(; lambda, beta)

    ## Solving
    niters = 10
    ub_step = 10
    solver = HiGHS.Optimizer

    ## The model
    function build_hydro(sp, t)
        @variable(sp, 0 <= vol <= 100, SDDP.State, initial_value = inivol)
        @variable(sp, 0 <= gh <= 60)
        @variable(sp, 0 <= spill <= 200)
        @variable(sp, 0 <= gt[i = 1:2] <= 15)
        @variable(sp, 0 <= def <= 75)
        @variable(sp, inflow)

        @constraint(sp, sum(gt) + gh + def == 75)
        @constraint(sp, vol.in + inflow - gh - spill == vol.out)

        if t == 1
            JuMP.fix(inflow, 0.0)
        else
            SDDP.parameterize(sp, inflows) do observed_inflow
                JuMP.fix(inflow, observed_inflow)
                return
            end
        end

        @stageobjective(sp, spill + 5 * gt[1] + 10 * gt[2] + 50 * def)
    end

    ## Lower bound via cuts
    println("Build primal outer model")
    pb = SDDP.LinearPolicyGraph(
        build_hydro;
        stages = nstages,
        sense = :Min,
        optimizer = solver,
        lower_bound = 0.0,
    )

    println("Solving primal outer model")
    SDDP.train(pb; iteration_limit = niters, risk_measure = rho)
    lb = SDDP.calculate_bound(pb; risk_measure = rho)
    println("       Risk-adjusted lower bound: ", round(lb; digits = 2))

    ## Monte-Carlo policy cost estimate; does not take risk_measure into account
    simulations = SDDP.simulate(pb, 500, [:vol, :gh, :spill, :gt, :def])
    objective_values =
        [sum(stage[:stage_objective] for stage in sim) for sim in simulations]
    μ = round(mean(objective_values); digits = 2)
    ci = round(1.96 * std(objective_values) / sqrt(500); digits = 2)

    println("Risk-neutral confidence interval: ", μ, " ± ", ci)

    ## Upper bound via inner approximation
    base_Lip = 50.0
    base_ub = 75.0 * base_Lip
    build_Lip(t) = base_Lip * (nstages - t)
    build_ub(t) = base_ub * (nstages - t)
    ibf = Inner.InnerBellmanFunction(
        build_Lip;
        upper_bound = build_ub,
        vertex_type = SDDP.SINGLE_CUT,
    )

    println("\nBuilding and Solving inner model for upper bounds")
    pb_inner, ub = Inner.inner_dp(
        build_hydro,
        pb;
        nstages,
        sense = :Min,
        optimizer = solver,
        lower_bound = 0.0,
        bellman_function = ibf,
        risk_measures = rho,
        print_level = 1,
    )

    @test lb <= ub
end

test_inner_hydro_1d()
