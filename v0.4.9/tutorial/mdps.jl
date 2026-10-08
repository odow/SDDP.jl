#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: Markov Decision Processes

# `SDDP.jl` can be used to solve a variety of Markov Decision processes. If the
# problem has continous state and control spaces, and the objective and
# transition function are convex, then SDDP.jl can find a globally optimal
# policy. In other cases, SDDP.jl will find a locally optimal policy.

# ## A simple example

# A simple demonstration of this is the example taken from page 98 of the book
# "Markov Decision Processes: Discrete stochastic Dynamic Programming", by
# Martin L. Putterman.

# The example, as described in Section 4.6.3 of the book, is to minimize a sum
# of squares of `N` non-negative variables, subject to a budget constraint that
# the variable values add up to `M`. Put mathematically, that is:
# ```math
# \begin{aligned}
# \min \;\; & \sum\limits_{i=1}^N x_i^2       \\
# s.t. \;\; & \sum\limits_{i=1}^N x_i = M     \\
#           & x_i \ge 0, \quad i \in 1,\ldots,N
# \end{aligned}
# ```
# The optimal objective value is ``M^2/N``, and the optimal solution is
# ``x_i = M / N``, which can be shown by induction.

# This can be reformulated as a Markov Decision Process by introducing a state
# variable, ``s``, which tracks the un-spent budget over ``N`` stages.

# ```math
# \begin{aligned}
# V_t(s) = \min \;\; & x^2 + V_{t+1}(s^\prime) \\
# s.t. \;\; & s^\prime = s - x \\
#           & x \le s \\
#           & x \ge 0 \\
#           & s \ge 0
# \end{aligned}
# ```
# and in the last stage ``V_N``, there is an additional constraint that
# ``s^\prime = 0``.

# The budget of ``M`` is computed by solving for ``V_1(M)``.

# !!! info
#     Since everything here is continuous and convex, SDDP.jl will find the
#     globally optimal policy.

# If the reformulation from the single problem into the recursive form of the
# Markov Decision Process is not obvious, consult Putterman's book.

# We can model and solve this problem using SDDP.jl as follows:

using SDDP
import Ipopt

M, N = 5, 3

model = SDDP.LinearPolicyGraph(
    stages = N,
    lower_bound = 0.0,
    optimizer = Ipopt.Optimizer,
) do subproblem, node
    @variable(subproblem, s >= 0, SDDP.State, initial_value = M)
    @variable(subproblem, x >= 0)
    @stageobjective(subproblem, x^2)
    @constraint(subproblem, x <= s.in)
    @constraint(subproblem, s.out == s.in - x)
    if node == N
        fix(s.out, 0.0; force = true)
    end
    return
end

SDDP.train(model; stopping_rules = [SDDP.BoundStalling(5, 1e-6)])

# Check that we got the theoretical optimum:

SDDP.calculate_bound(model), M^2 / N

# And check that we cound the theoretical value for each ``x_i``:

simulations = SDDP.simulate(model, 1, [:x])
for data in simulations[1]
    println("x_$(data[:node_index]) = $(data[:x])")
end

# Close enough! We don't get exactly 5/3 because of numerical tolerances within
# our choice of optimization solver (in this case, Ipopt).

# ## A more complicated policy

# SDDP.jl is also capable of finding policies for other types of Markov Decision
# Processes. A classic example of a Markov Decision Process is the problem of
# finding a path through a maze.

# Here's one example of a maze. Try changing the parameters to explore different
# mazes:

M, N = 3, 4
initial_square = (1, 1)
reward, illegal_squares, penalties = (3, 4), [(2, 2)], [(3, 1), (2, 4)]
path = fill("⋅", M, N)
path[initial_square...] = "1"
for (k, v) in (illegal_squares => "▩", penalties => "†", [reward] => "*")
    for (i, j) in k
        path[i, j] = v
    end
end
print(join([join(path[i, :], ' ') for i in 1:size(path, 1)], '\n'))

# Our goal is to get from square `1` to square `*`. If we step on a `†`, we
# incur a penalty of `1`. Squares with `▩` are blocked; we cannot move there.

# There are a variety of ways that we can solve this problem. We're going to
# solve it using a stationary binary stochastic programming formulation.

# Our state variable will be a matrix of binary variables ``x_{i,j}``, where
# each element is ``1`` if the agent is in the square and ``0`` otherwise. In
# each period, we incur a reward of ``1`` if we are in the `reward` square and a
# penalty of ``-1`` if we are in a `penalties` square. We cannot move to the
# `illegal_squares`, so those ``x_{i,j} = 0``. Feasibility between moves is
# modelled by constraints of the form:
# ```math
# x^\prime_{i,j} \le \sum\limits_{(a,b)\in P} x_{a,b}
# ```
# where ``P`` is the set of squares from which it is valid to move from `(a, b)`
# to `(i, j)`.

# Because we are looking for a stationary policy, we need a unicyclic graph with
# a discount factor:

discount_factor = 0.9
graph = SDDP.UnicyclicGraph(discount_factor)

# Then we can formulate our full model:

import HiGHS

model = SDDP.PolicyGraph(
    graph;
    sense = :Max,
    upper_bound = 1 / (1 - discount_factor),
    optimizer = HiGHS.Optimizer,
) do sp, _
    ## Our state is a binary variable for each square
    @variable(
        sp,
        x[i = 1:M, j = 1:N],
        Bin,
        SDDP.State,
        initial_value = (i, j) == initial_square,
    )
    ## Can only be in one square at a time
    @constraint(sp, sum(x[i, j].out for i in 1:M, j in 1:N) == 1)
    ## Incur rewards and penalties
    @stageobjective(
        sp,
        x[reward...].out - sum(x[i, j].out for (i, j) in penalties)
    )
    ## Some squares are illegal
    @constraint(sp, [(i, j) in illegal_squares], x[i, j].out <= 0)
    ## Constraints on valid moves
    for i in 1:M, j in 1:N
        moves = [(i - 1, j), (i + 1, j), (i, j), (i, j + 1), (i, j - 1)]
        filter!(v -> 1 <= v[1] <= M && 1 <= v[2] <= N, moves)
        @constraint(sp, x[i, j].out <= sum(x[a, b].in for (a, b) in moves))
    end
    return
end

# The upper bound is obtained by assuming that we reach the reward square in one
# move and stay there.

# !!! warning
#     Since there are discrete decisions here, SDDP.jl is not guaranteed to find
#     the gobally optimal policy.

SDDP.train(model; stopping_rules = [SDDP.BoundStalling(5, 1e-6)])

# Simulating a cyclic policy graph requires an explicit `sampling_scheme` that
# does not terminate early based on the cycle probability:

simulations = SDDP.simulate(
    model,
    1,
    [:x];
    sampling_scheme = SDDP.InSampleMonteCarlo(
        max_depth = 5,
        terminate_on_dummy_leaf = false,
    ),
);

# Fill in the `path` with the time-step in which we visit the square:

for (t, data) in enumerate(simulations[1]), i in 1:M, j in 1:N
    if data[:x][i, j].in > 0.5
        path[i, j] = "$t"
    end
end

print(join([join(path[i, :], ' ') for i in 1:size(path, 1)], '\n'))

# !!! tip
#     This formulation will likely struggle as the number of cells in the maze
#     increases. Can you think of an equivalent formulation that uses fewer
#     state variables?
