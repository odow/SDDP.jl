#  Copyright 2017-21, Oscar Dowson.                                     #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Biobjective hydro-thermal

using SDDP, GLPK, Statistics, Test

function biobjective_example()
    model = SDDP.LinearPolicyGraph(
        stages = 3,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do subproblem, _
        @variable(subproblem, 0 <= v <= 200, SDDP.State, initial_value = 50)
        @variables(subproblem, begin
            0 <= g[i = 1:2] <= 100
            0 <= u <= 150
            s >= 0
            shortage_cost >= 0
        end)
        @expressions(subproblem, begin
            objective_1, g[1] + 10 * g[2]
            objective_2, shortage_cost
        end)
        @constraints(subproblem, begin
                inflow_constraint, v.out == v.in - u - s
                g[1] + g[2] + u == 150
                shortage_cost >= 40 - v.out
                shortage_cost >= 60 - 2 * v.out
                shortage_cost >= 80 - 4 * v.out
            end)
        ## You must call this for a biobjective problem!
        SDDP.initialize_biobjective_subproblem(subproblem)
        SDDP.parameterize(subproblem, 0.0:5:50.0) do ω
            JuMP.set_normalized_rhs(inflow_constraint, ω)
            ## You must call `set_biobjective_functions` from within
            ## `SDDP.parameterize`.
            return SDDP.set_biobjective_functions(
                subproblem,
                objective_1,
                objective_2,
            )
        end
    end
    pareto_weights =
        SDDP.train_biobjective(model, solution_limit = 10, iteration_limit = 10)
    solutions = [(k, v) for (k, v) in pareto_weights]
    sort!(solutions; by = x -> x[1])
    @test length(solutions) == 10
    ## Test for convexity! The gradient must be decreasing as we move from left
    ## to right.
    gradient(a, b) = (b[2] - a[2]) / (b[1] - a[1])
    grad = Inf
    for i in 1:9
        new_grad = gradient(solutions[i], solutions[i+1])
        @test new_grad < grad
        grad = new_grad
    end
    return
end

biobjective_example()
