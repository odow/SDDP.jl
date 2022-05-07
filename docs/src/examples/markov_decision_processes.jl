#  Copyright (c) 2017-22, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Markov Decision Processes

# This example is taken from page 98 of the book
# [Markov Decision Processes: Discrete stochastic Dynamic Programming](https://onlinelibrary.wiley.com/doi/book/10.1002/9780470316887),
# by Martin L. Putterman.

# The example is described in section 4.6.3 of the book as:
# ```math
# \begin{aligned}
# \min \;     & \sum\limits_{i = 1}^N x_i^2            \\
# \text{s.t.} & \sum\limits_{i = 1}^N x_i = M          \\
#             & x_i \ge 0 \quad \forall i = 1,\ldots,N.
# \end{aligned}
# ```

# The optimal objective value is ``\frac{M^2}{N}``, and the optimal primal
# solution is ``x_i = \frac{M}{N}, i=1,\ldots,N``, which can be shown by
# induction.

# In this tutorial, we demonstrate two ways to solve this problem. First as a
# single nonlinear optimization problem, and second using SDDP.jl.

using SDDP
using JuMP
using Test

import Ipopt

# ## Nonlinear program

function solve_as_nonlinear(; N::Int, M::Int)
    model = Model(Ipopt.Optimizer)
    @variable(model, x[1:N] >= 0)
    @constraint(model, sum(x[i] for i in 1:N) == M)
    @objective(model, Min, sum(x[i]^2 for i in 1:N))
    optimize!(model)
    return value.(x)
end

x = solve_as_nonlinear(; N = 3, M = 5)

#-

@test all(isapprox.(x, M / N; atol = 1e-2))

# ## Using SDDP.jl

function solve_as_sddp(; N::Int, M::Int)
    model = SDDP.LinearPolicyGraph(
        stages = N,
        lower_bound = 0.0,
        optimizer = Ipopt.Optimizer,
    ) do subproblem, node
        @variable(subproblem, x >= 0)
        @stageobjective(subproblem, x^2)
        S = (node == N) ? 0 : M
        @variable(subproblem, 0 <= s <= S, SDDP.State, initial_value = M)
        @constraint(subproblem, s.out == s.in - x)
    end
    SDDP.train(model, iteration_limit = 30)
    ## Now we need to extract the primal solution. Simulate the policy:
    simulations = SDDP.simulate(model, 1, [:x])
    ## And extract `x` from the first simulation:
    return [simulation[:x] for simulation in simulations[1]]
end

x = solve_as_sddp(; N = 3, M = 5)

#-

@test all(isapprox.(x, M / N; atol = 1e-2))
