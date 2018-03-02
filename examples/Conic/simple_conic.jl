#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Gurobi, Base.Test

"""
A trivial two-stage stochastic programming problem to work out the
least-squares estimator for the set of `data`. The first stage chooses `x`,
and the second-stage evaluates the decision.
"""
function least_squares(data::Vector{Float64})
    m = SDDPModel(
        stages = 2,
        sense = :Min,
        objective_bound = 0.0,
        solver=GurobiSolver(OutputFlag=0)
            ) do sp, t

        @state(sp, x′>=0, x==0)
        if t == 2
            @stageobjective(sp, ω=data, (x - ω)^2)
        else
            @stageobjective(sp, 0.0)
        end
    end
end

srand(1234)
data = rand(50)
m = least_squares(data)
status = solve(m, max_iterations=20, print_level=0)
@test status == :max_iterations
@test isapprox(getbound(m), mean((data - mean(data)).^2), atol=1e-6)
