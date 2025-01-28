#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: two-stage newsvendor

# The purpose of this tutorial is to demonstrate how to model and solve a
# two-stage stochastic program.

# It is based on the [Two stage stochastic programs](https://jump.dev/JuMP.jl/dev/tutorials/applications/two_stage_stochastic/)
# tutorial in JuMP.

# This tutorial uses the following packages

using JuMP
using SDDP
import Distributions
import ForwardDiff
import HiGHS
import Plots
import StatsPlots
import Statistics

# ## Background

# The data for this problem is:

D = Distributions.TriangularDist(150.0, 250.0, 200.0)
N = 100
d = sort!(rand(D, N));
Ω = 1:N
P = fill(1 / N, N);
StatsPlots.histogram(d; bins = 20, label = "", xlabel = "Demand")

# ## Kelley's cutting plane algorithm

# Kelley's cutting plane algorithm is an iterative method for maximizing concave
# functions. Given a concave function $f(x)$, Kelley's constructs an
# outer-approximation of the function at the minimum by a set of first-order
# Taylor series approximations (called **cuts**) constructed at a set of points
# $k = 1,\ldots,K$:
# ```math
# \begin{aligned}
# f^K = \max\limits_{\theta \in \mathbb{R}, x \in \mathbb{R}^N} \;\; & \theta\\
# & \theta \le f(x_k) + \nabla f(x_k)^\top (x - x_k),\quad k=1,\ldots,K\\
# & \theta \le M,
# \end{aligned}
# ```
# where $M$ is a sufficiently large number that is an upper bound for $f$ over
# the domain of $x$.

# Kelley's cutting plane algorithm is a structured way of choosing points $x_k$
# to visit, so that as more cuts are added:
# ```math
# \lim_{K \rightarrow \infty} f^K = \max\limits_{x \in \mathbb{R}^N} f(x)
# ```
# However, before we introduce the algorithm, we need to introduce some bounds.

# ### Bounds

# By convexity, $f(x) \le f^K$ for all $x$. Thus, if $x^*$ is a maximizer of
# $f$, then at any point in time we can construct an upper bound for $f(x^*)$ by
# solving $f^K$.

# Moreover, we can use the primal solutions $x_k^*$ returned by solving $f^k$ to
# evaluate $f(x_k^*)$ to generate a lower bound.

# Therefore, $\max\limits_{k=1,\ldots,K} f(x_k^*) \le f(x^*) \le f^K$.

# When the lower bound is sufficiently close to the upper bound, we can
# terminate the algorithm and declare that we have found an solution that is
# close to optimal.

# ### Implementation

# Here is pseudo-code fo the Kelley algorithm:

# 1. Take as input a convex function $f(x)$ and a iteration limit $K_{max}$.
#    Set $K = 1$, and initialize $f^{K-1}$. Set $lb = -\infty$ and $ub = \infty$.
# 2. Solve $f^{K-1}$ to obtain a candidate solution $x_{K}$.
# 3. Update $ub = f^{K-1}$ and $lb = \max\{lb, f(x_{K})\}$.
# 4. Add a cut $\theta \ge f(x_{K}) + \nabla f\left(x_{K}\right)^\top (x - x_{K})$ to form $f^{K}$.
# 5. Increment $K$.
# 6. If $K > K_{max}$ or $|ub - lb| < \epsilon$, STOP, otherwise, go to step 2.

# And here's a complete implementation:

function kelleys_cutting_plane(
    ## The function to be minimized.
    f::Function,
    ## The gradient of `f`. By default, we use automatic differentiation to
    ## compute the gradient of f so the user doesn't have to!
    ∇f::Function = x -> ForwardDiff.gradient(f, x);
    ## The number of arguments to `f`.
    input_dimension::Int,
    ## An upper bound for the function `f` over its domain.
    upper_bound::Float64,
    ## The number of iterations to run Kelley's algorithm for before stopping.
    iteration_limit::Int,
    ## The absolute tolerance ϵ to use for convergence.
    tolerance::Float64 = 1e-6,
)
    ## Step (1):
    K = 1
    model = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(model)
    JuMP.@variable(model, θ <= upper_bound)
    JuMP.@variable(model, x[1:input_dimension])
    JuMP.@objective(model, Max, θ)
    x_k = fill(NaN, input_dimension)
    lower_bound, upper_bound = -Inf, Inf
    while true
        ## Step (2):
        JuMP.optimize!(model)
        x_k .= JuMP.value.(x)
        ## Step (3):
        upper_bound = JuMP.objective_value(model)
        lower_bound = min(upper_bound, f(x_k))
        println("K = $K : $(lower_bound) <= f(x*) <= $(upper_bound)")
        ## Step (4):
        JuMP.@constraint(model, θ <= f(x_k) + ∇f(x_k)' * (x .- x_k))
        ## Step (5):
        K = K + 1
        ## Step (6):
        if K > iteration_limit
            println("-- Termination status: iteration limit --")
            break
        elseif abs(upper_bound - lower_bound) < tolerance
            println("-- Termination status: converged --")
            break
        end
    end
    println("Found solution: x_K = ", x_k)
    return
end

# Let's run our algorithm to see what happens:

kelleys_cutting_plane(;
    input_dimension = 2,
    upper_bound = 10.0,
    iteration_limit = 20,
) do x
    return -(x[1] - 1)^2 + -(x[2] + 2)^2 + 1.0
end

# ## L-Shaped theory

# The L-Shaped method is a way of solving two-stage stochastic programs by
# Benders' decomposition. It takes the problem:

# ```math
# \begin{aligned}
# V = \max\limits_{x,y_\omega} \;\; & -2x + \mathbb{E}_\omega[5y_\omega - 0.1(x - y_\omega)] \\
#   & y_\omega \le x              & \quad \forall \omega \in \Omega \\
#   & 0 \le y_\omega \le d_\omega & \quad \forall \omega \in \Omega \\
#   & x \ge 0.
# \end{aligned}
# ```

# and decomposes it into a second-stage problem:

# ```math
# \begin{aligned}
# V_2(\bar{x}, d_\omega) = \max\limits_{x,x^\prime,y_\omega} \;\; & 5y_\omega - x^\prime \\
#   & y_\omega \le x \\
#   & x^\prime = x - y_\omega \\
#   & 0 \le y_\omega \le d_\omega \\
#   & x = \bar{x} & [\lambda]
# \end{aligned}
# ```

# and a first-stage problem:

# ```math
# \begin{aligned}
# V = \max\limits_{x,\theta} \;\; & -2x + \theta \\
#   & \theta \le \mathbb{E}_\omega[V_2(x, \omega)] \\
#   & x \ge 0
# \end{aligned}
# ```

# Then, because $V_2$ is convex with respect to $\bar{x}$ for fixed $\omega$,
# we can use a set of feasible points $\{x^k\}$ construct an outer approximation:
# ```math
# \begin{aligned}
# V^K = \max\limits_{x,\theta} \;\; & -2x + \theta \\
#   & \theta \le \mathbb{E}_\omega[V_2(x^k, \omega) + \nabla V_2(x^k, \omega)^\top(x - x^k)] & \quad k = 1,\ldots,K\\
#   & x \ge 0 \\
#   & \theta \le M
# \end{aligned}
# ```
# where $M$ is an upper bound on possible values of $V_2$ so that the problem
# has a bounded solution.

# It is also useful to see that because $\bar{x}$ appears only on the right-hand
# side of a linear program, $\nabla V_2(x^k, \omega) = \lambda^k$.

# Ignoring how we choose $x^k$ for now, we can construct a lower and upper bound
# on the optimal solution:

# $$-2x^K + \mathbb{E}_\omega[V_2(x^K, \omega)] = \underbar{V} \le V \le \overline{V} = V^K$$

# Thus, we need some way of cleverly choosing a sequence of $x^k$ so that the
# lower bound converges to the upper bound.

# 1. Start with $K=1$
# 2. Solve $V^{K-1}$ to get $x^K$
# 3. Set $\overline{V} = V^k$
# 4. Solve $V_2(x^K, \omega)$ for all $\omega$ and store the optimal objective
#    value and dual solution $\lambda^K$
# 5. Set $\underbar{V} = -2x^K + \mathbb{E}_\omega[V_2(x^k, \omega)]$
# 6. If $\underbar{V} \approx \overline{V}$, STOP
# 7. Add new constraint $\theta \le \mathbb{E}_\omega[V_2(x^K, \omega) +\lambda^K (x - x^K)]$
# 8. Increment $K$, GOTO 2

# The next section implements this algorithm in Julia.

# ## L-Shaped implementation

# Here's a function to compute the second-stage problem;

function solve_second_stage(x̅, d_ω)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x_in)
    @variable(model, x_out >= 0)
    fix(x_in, x̅)
    @variable(model, 0 <= u_sell <= d_ω)
    @constraint(model, x_out == x_in - u_sell)
    @constraint(model, u_sell <= x_in)
    @objective(model, Max, 5 * u_sell - 0.1 * x_out)
    optimize!(model)
    return (
        V = objective_value(model),
        λ = reduced_cost(x_in),
        x = value(x_out),
        u = value(u_sell),
    )
end

solve_second_stage(200, 170)

# Here's the first-stage subproblem:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x_in == 0)
@variable(model, x_out >= 0)
@variable(model, u_make >= 0)
@constraint(model, x_out == x_in + u_make)
M = 5 * maximum(d)
@variable(model, θ <= M)
@objective(model, Max, -2 * u_make + θ)

# Importantly, to ensure we have a bounded solution, we need to add an upper
# bound to the variable `θ`.

kIterationLimit = 100
for k in 1:kIterationLimit
    println("Solving iteration k = $k")
    ## Step 2
    optimize!(model)
    xᵏ = value(x_out)
    println("  xᵏ = $xᵏ")
    ## Step 3
    ub = objective_value(model)
    println("  V̅ = $ub")
    ## Step 4
    ret = [solve_second_stage(xᵏ, d[ω]) for ω in Ω]
    ## Step 5
    lb = value(-2 * u_make) + sum(p * r.V for (p, r) in zip(P, ret))
    println("  V̲ = $lb")
    ## Step 6
    if ub - lb < 1e-6
        println("Terminating with near-optimal solution")
        break
    end
    ## Step 7
    c = @constraint(
        model,
        θ <= sum(p * (r.V + r.λ * (x_out - xᵏ)) for (p, r) in zip(P, ret)),
    )
    println("  Added cut: $c")
end

# To get the first-stage solution, we do:

optimize!(model)
xᵏ = value(x_out)

# To compute a second-stage solution, we do:

solve_second_stage(xᵏ, 170.0)

# ## Policy Graph

# Now let's see how we can formulate and train a policy for the two-stage
# newsvendor problem using `SDDP.jl`. Under the hood, `SDDP.jl` implements the
# exact algorithm that we just wrote by hand.

model = SDDP.LinearPolicyGraph(;
    stages = 2,
    sense = :Max,
    upper_bound = 5 * maximum(d),  # The `M` in θ <= M
    optimizer = HiGHS.Optimizer,
) do subproblem::JuMP.Model, stage::Int
    @variable(subproblem, x >= 0, SDDP.State, initial_value = 0)
    if stage == 1
        @variable(subproblem, u_make >= 0)
        @constraint(subproblem, x.out == x.in + u_make)
        @stageobjective(subproblem, -2 * u_make)
    else
        @variable(subproblem, u_sell >= 0)
        @constraint(subproblem, u_sell <= x.in)
        @constraint(subproblem, x.out == x.in - u_sell)
        SDDP.parameterize(subproblem, d, P) do ω
            set_upper_bound(u_sell, ω)
            return
        end
        @stageobjective(subproblem, 5 * u_sell - 0.1 * x.out)
    end
    return
end

SDDP.train(model; log_every_iteration = true)

# One way to query the optimal policy is with [`SDDP.DecisionRule`](@ref):

first_stage_rule = SDDP.DecisionRule(model; node = 1)

#-

solution_1 = SDDP.evaluate(first_stage_rule; incoming_state = Dict(:x => 0.0))

# Here's the second stage:

second_stage_rule = SDDP.DecisionRule(model; node = 2)
solution = SDDP.evaluate(
    second_stage_rule;
    incoming_state = Dict(:x => solution_1.outgoing_state[:x]),
    noise = 170.0,  # A value of d[ω], can be out-of-sample.
    controls_to_record = [:u_sell],
)

# ## Simulation

# Querying the decision rules is tedious. It's often more useful to simulate the
# policy:

simulations = SDDP.simulate(
    model,
    10,  #= number of replications =#
    [:x, :u_sell, :u_make];  #= variables to record =#
    skip_undefined_variables = true,
);

# `simulations` is a vector with 10 elements

length(simulations)

# and each element is a vector with two elements (one for each stage)

length(simulations[1])

# The first stage contains:

simulations[1][1]

# The second stage contains:

simulations[1][2]

# We can compute aggregated statistics across the simulations:

objectives = map(simulations) do simulation
    return sum(data[:stage_objective] for data in simulation)
end
μ, t = SDDP.confidence_interval(objectives)
println("Simulation ci : $μ ± $t")

# ## Risk aversion revisited

# SDDP.jl contains a number of risk measures. One example is:

0.5 * SDDP.Expectation() + 0.5 * SDDP.WorstCase()

# You can construct a risk-averse policy by passing a risk measure to the
# `risk_measure` keyword argument of [`SDDP.train`](@ref).

# We can explore how the optimal decision changes with risk by creating a
# function:

function solve_newsvendor(risk_measure::SDDP.AbstractRiskMeasure)
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        sense = :Max,
        upper_bound = 5 * maximum(d),
        optimizer = HiGHS.Optimizer,
    ) do subproblem, node
        @variable(subproblem, x >= 0, SDDP.State, initial_value = 0)
        if node == 1
            @stageobjective(subproblem, -2 * x.out)
        else
            @variable(subproblem, u_sell >= 0)
            @constraint(subproblem, u_sell <= x.in)
            @constraint(subproblem, x.out == x.in - u_sell)
            SDDP.parameterize(subproblem, d, P) do ω
                set_upper_bound(u_sell, ω)
                return
            end
            @stageobjective(subproblem, 5 * u_sell - 0.1 * x.out)
        end
        return
    end
    SDDP.train(model; risk_measure = risk_measure, print_level = 0)
    first_stage_rule = SDDP.DecisionRule(model; node = 1)
    solution = SDDP.evaluate(first_stage_rule; incoming_state = Dict(:x => 0.0))
    return solution.outgoing_state[:x]
end

# Now we can see how many units a decision maker would order using `CVaR`:

solve_newsvendor(SDDP.CVaR(0.4))

# as well as a decision-maker who cares only about the worst-case outcome:

solve_newsvendor(SDDP.WorstCase())

# In general, the decision-maker will be somewhere between the two extremes.
# The [`SDDP.Entropic`](@ref) risk measure is a risk measure that has a single
# parameter that lets us explore the space of policies between the two extremes.
# When the parameter is small, the measure acts like [`SDDP.Expectation`](@ref),
# and when it is large, it acts like [`SDDP.WorstCase`](@ref).

# Here is what we get if we solve our problem multiple times for different
# values of the risk aversion parameter ``\gamma``:

Γ = [10^i for i in -4:0.5:1]
buy = [solve_newsvendor(SDDP.Entropic(γ)) for γ in Γ]
Plots.plot(
    Γ,
    buy;
    xaxis = :log,
    xlabel = "Risk aversion parameter γ",
    ylabel = "Number of pies to make",
    legend = false,
)

# ## Things to try

# There are a number of things you can try next:

#  * Experiment with different buy and sales prices
#  * Experiment with different distributions of demand
#  * Explore how the optimal policy changes if you use a different risk measure
#  * What happens if you can only buy and sell integer numbers of newspapers?
#    Try this by adding `Int` to the variable definitions:
#    `@variable(subproblem, buy >= 0, Int)`
#  * What happens if you use a different upper bound? Try an invalid one like
#    `-100`, and a very large one like `1e12`.
