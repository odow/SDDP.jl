#  Copyright 2017-20, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#==
Modified version of the Asset Management problem taken from

    J. R. Birge,  F. Louveaux,  Introduction to Stochastic Programming,
    Springer Series in Operations Research and Financial Engineering,
    Springer New York, New York, NY, 2011
==#

using SDDP, GLPK, Test

function asset_management_stagewise(; cut_type)
    ws = [1.25, 1.06]
    wb = [1.14, 1.12]
    Phi = [-1, 5]
    Psi = [0.02, 0.0]

    model = SDDP.MarkovianPolicyGraph(
        sense = :Max,
        transition_matrices = Array{Float64,2}[
            [1.0]',
            [0.5 0.5],
            [0.5 0.5; 0.5 0.5],
            [0.5 0.5; 0.5 0.5],
        ],
        bellman_function = SDDP.BellmanFunction(upper_bound = 1000.0, cut_type = cut_type),
        optimizer = GLPK.Optimizer,
    ) do subproblem, node
        t, i = node
        @variable(subproblem, xs >= 0, SDDP.State, initial_value = 0)
        @variable(subproblem, xb >= 0, SDDP.State, initial_value = 0)
        if t == 1
            @constraint(subproblem, xs.out + xb.out == 55 + xs.in + xb.in)
            @stageobjective(subproblem, 0)
        elseif t == 2 || t == 3
            @variable(subproblem, phi)
            @constraint(subproblem, ws[i] * xs.in + wb[i] * xb.in + phi == xs.out + xb.out)
            SDDP.parameterize(subproblem, [1, 2], [0.6, 0.4]) do ω
                JuMP.fix(phi, Phi[ω])
                @stageobjective(subproblem, Psi[ω] * xs.out)
            end
        else
            @variable(subproblem, u >= 0)
            @variable(subproblem, v >= 0)
            @constraint(subproblem, ws[i] * xs.in + wb[i] * xb.in + u - v == 80)
            @stageobjective(subproblem, -4u + v)
        end
    end
    SDDP.train(
        model;
        iteration_limit = 100,
        print_level = 0,
        risk_measure = (node) -> begin
            if node[1] != 3
                SDDP.Expectation()
            else
                SDDP.EAVaR(lambda = 0.5, beta = 0.5)
            end
        end,
    )
    @test SDDP.calculate_bound(model) ≈ 1.278 atol = 1e-3
end

asset_management_stagewise(cut_type = SDDP.SINGLE_CUT)
asset_management_stagewise(cut_type = SDDP.MULTI_CUT)
