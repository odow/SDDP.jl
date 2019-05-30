#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#==
    This example comes from
        https://github.com/blegat/StochasticDualDynamicProgramming.jl/blob/fe5ef82db6befd7c8f11c023a639098ecb85737d/test/prob5.2_3stages.jl
==#

using SDDP, GLPK, Test

function test_prob52_3stages()
    model = SDDP.LinearPolicyGraph(
                stages = 3,
                lower_bound = 0.0,
                optimizer = with_optimizer(GLPK.Optimizer)) do sp, t
        n = 4
        m = 3
        ic = [16, 5, 32, 2]
        C = [25, 80, 6.5, 160]
        T = [8760, 7000, 1500] / 8760
        D2 = [diff([0, 3919, 7329, 10315])  diff([0, 7086, 9004, 11169])]
        p2 = [0.9, 0.1]
        @variable(sp, x[i=1:n] >= 0, SDDP.State, initial_value = 0.0)
        @variables(sp, begin
            y[1:n, 1:m] >= 0
            v[1:n]      >= 0
            penalty     >= 0
            ξ[j=1:m]
        end)
        @constraints(sp, begin
            [i=1:n], x[i].out == x[i].in + v[i]
            [i=1:n], sum(y[i, :]) <= x[i].in
            [j=1:m], sum(y[:, j]) + penalty >= ξ[j]
        end)
        @stageobjective(sp, ic'v + C' * y * T + 1e5 * penalty)
        if t != 1 # no uncertainty in first stage
            SDDP.parameterize(sp, 1:size(D2, 2), p2) do ω
                for j in 1:m
                    JuMP.fix(ξ[j], D2[j, ω])
                end
            end
        end
        if t == 3
            @constraint(sp, sum(v) == 0)
        end
    end

    det = SDDP.deterministic_equivalent(model, with_optimizer(GLPK.Optimizer))
    JuMP.optimize!(det)
    @test JuMP.objective_value(det) ≈ 406712.49 atol = 0.1

    SDDP.train(model, iteration_limit = 30, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 406712.49 atol = 0.1
end

test_prob52_3stages()
