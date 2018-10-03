#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==
    An implementation of the Production Management example from FAST
    https://github.com/leopoldcambier/FAST/blob/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/production management multiple stages/
==#

using Kokako, GLPK, Test

function fast_production_management()
    DEMAND = [2, 10]
    H = 3
    N = 2
    C = [0.2, 0.7]
    S = 2 .+ [0.33, 0.54]
    model = Kokako.PolicyGraph(Kokako.LinearGraph(H),
                bellman_function = Kokako.AverageCut(lower_bound=-50.0),
                optimizer = with_optimizer(GLPK.Optimizer)
                        ) do sp, t
        @variables(sp, begin
            x0[i=1:N]
            x[i=1:N] >= 0
            s[i=1:N] >= 0
            d
        end)
        Kokako.add_state_variable.(sp, x0, x, 0.0)
        @constraints(sp, begin
            s .<= x0
            sum(s) <= d
        end)
        Kokako.parameterize(sp, t==1 ? [0] : DEMAND) do ω
            JuMP.fix(d, ω)
        end
        @stageobjective(sp, C'x - S's)
    end
    Kokako.train(model, iteration_limit = 10, print_level = 0)
    @test Kokako.calculate_bound(model) ≈ -23.96 atol=1e-2
end

fast_production_management()
