#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using GLPK
using SDDP
using Test

function _create_model(minimization::Bool)
    DEMAND = [2, 10]
    H = 3
    N = 2
    C = [0.2, 0.7]
    S = 2 .+ [0.33, 0.54]
    sense = minimization ? :Min : :Max
    model = SDDP.LinearPolicyGraph(
        stages = H,
        sense = sense,
        lower_bound = -50.0,
        upper_bound = 50.0,
    ) do sp, t
        @variable(sp, x[1:N] >= 0, SDDP.State, initial_value = 0.0)
        @variables(sp, begin
            s[i = 1:N] >= 0
            d
        end)
        @constraints(sp, begin
            [i = 1:N], s[i] <= x[i].in
            c, sum(s) <= d + 1
        end)
        SDDP.parameterize(sp, t == 1 ? [1] : 1:length(DEMAND)) do ω
            JuMP.fix(d, DEMAND[ω])
            set_upper_bound(s[1], 0.1 * ω)
            set_lower_bound(x[1].out, ω)
            set_normalized_rhs(c, ω)
            sgn = minimization ? 1.0 : -1.0
            @stageobjective(
                sp, sgn * (sum(C[i] * x[i].out for i = 1:N) - S[ω] * s[ω] + ω)
            )
        end
    end
    return model
end

@testset "Roundtrips" begin
    @testset "Min: Read and write to file" begin
        base_model = _create_model(true)
        set_optimizer(base_model, GLPK.Optimizer)
        SDDP.train(base_model; iteration_limit = 50, print_level = 0)

        model = _create_model(true)
        SDDP.write_to_file(model, "experimental.sof.json")
        set_optimizer(model, GLPK.Optimizer)
        SDDP.train(model; iteration_limit = 50, print_level = 0)

        new_model = SDDP.read_from_file("experimental.sof.json")
        set_optimizer(new_model, GLPK.Optimizer)
        SDDP.train(new_model; iteration_limit = 50, print_level = 0)

        @test isapprox(
            SDDP.calculate_bound(base_model),
            SDDP.calculate_bound(model);
            atol = 1e-6
        )

        @test isapprox(
            SDDP.calculate_bound(base_model),
            SDDP.calculate_bound(new_model);
            atol = 1e-6
        )
    end

    @testset "Max: Read and write to file" begin
        base_model = _create_model(false)
        set_optimizer(base_model, GLPK.Optimizer)
        SDDP.train(base_model; iteration_limit = 50, print_level = 0)

        model = _create_model(false)
        SDDP.write_to_file(model, "experimental.sof.json")
        set_optimizer(model, GLPK.Optimizer)
        SDDP.train(model; iteration_limit = 50, print_level = 0)

        new_model = SDDP.read_from_file("experimental.sof.json")
        set_optimizer(new_model, GLPK.Optimizer)
        SDDP.train(new_model; iteration_limit = 50, print_level = 0)

        @test isapprox(
            SDDP.calculate_bound(base_model),
            SDDP.calculate_bound(model);
            atol = 1e-6
        )

        @test isapprox(
            SDDP.calculate_bound(base_model),
            SDDP.calculate_bound(new_model);
            atol = 1e-6
        )
    end
    rm("experimental.sof.json")
end
