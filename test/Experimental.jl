#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP
using GLPK

function FAST_production_scheduling()
    DEMAND = [2, 10]
    H = 3
    N = 2
    C = [0.2, 0.7]
    S = 2 .+ [0.33, 0.54]
    model = SDDP.LinearPolicyGraph(
        stages = H,
        lower_bound = -50.0,
    ) do sp, t
        @variable(sp, x[1:N] >= 0, SDDP.State, initial_value = 0.0)
        @variables(sp, begin
            s[i = 1:N] >= 0
            d
        end)
        @constraints(sp, begin
            [i = 1:N], s[i] <= x[i].in
            sum(s) <= d
        end)
        SDDP.parameterize(sp, t == 1 ? [0] : DEMAND) do ω
            JuMP.fix(d, ω)
        end
        @stageobjective(sp, sum(C[i] * x[i].out for i = 1:N) - S's)
    end
    return model
end

@testset "Read and write to file" begin
    model = FAST_production_scheduling()
    SDDP.Experimental.write_to_file(model, "experimental.sof.json")
    new_model = SDDP.Experimental.read_from_file("experimental.sof.json")
    set_optimizer(new_model, GLPK.Optimizer)
    SDDP.train(new_model; iteration_limit = 50, print_level = 0)
    @test SDDP.calculate_bound(new_model) ≈ -23.96 atol = 1e-2
    rm("experimental.sof.json")
end
