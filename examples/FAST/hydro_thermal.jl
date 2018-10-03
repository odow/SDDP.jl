#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==
    An implementation of the Hydro-thermal example from FAST
    https://github.com/leopoldcambier/FAST/tree/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/hydro%20thermal
==#

using Kokako, GLPK, Test

function fast_hydro_thermal()
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                bellman_function = Kokako.AverageCut(lower_bound=0.0),
                optimizer = with_optimizer(GLPK.Optimizer)
                        ) do sp, t
        # @state(sp, 0 <= x <= 8, x0==0)
        @variables(sp, begin
            x0
            0 <= x <= 8
            y >= 0
            p >= 0
            ξ
        end)
        Kokako.add_state_variable(sp, x0, x, 0.0)
        @constraints(sp, begin
            p + y >= 6
            x <= x0 - y + ξ
        end)
        RAINFALL = (t == 1 ? [6] : [2, 10])
        Kokako.parameterize(sp, RAINFALL) do ω
            JuMP.fix(ξ, ω)
        end
        @stageobjective(sp, 5 * p)
    end

    Kokako.train(model, iteration_limit = 10, print_level = 0)

    @test Kokako.calculate_bound(model) == 10
end

fast_hydro_thermal()
