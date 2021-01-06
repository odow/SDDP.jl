# # Expert I: Pedagogical SDDP

# In this tutorial we walk through a simplified implementation of stochastic
# dual dynamic programming to explain the key concepts.

# !!! note
#     If you haven't already, go read [Basic I: first steps](@ref), since it
#     introduces much necessary background theory, and
#     [Basic II: adding uncertainty](@ref), since it defines the same example we
#     will solve in our implementation.

# This tutorial uses the following packages. For clarity, we call
# `import PackageName` so that we must prefix `PackageName.` to all functions
# and structs provided by that package. Everything not prefixed is either part
# of base Julia, or we wrote it.

import ForwardDiff
import GLPK
import JuMP
import Statistics

# ## Theory

# ### Preliminaries: Kelley's cutting plane algorithm

# Kelley's cutting plane algorithm is an iterative method for minimizing convex
# functions. Given a convex function $f(x)$, Kelley's constructs an
# under-approximation of the function at the minimum by a set of first-order
# Taylor series approximations (called **cuts**) constructed at a set of $K$
# points $k = 1,\ldots,K$:
# ```math
# \begin{aligned}
# f^K = \min\limits_{\theta, x} \;\; & \theta\\
# & \theta \ge f(x_k) + \frac{df}{dx}\left(x_k\right) \cdot (x - x_k),\quad k=1,\ldots,K\\
# & \theta \ge M,
# \end{aligned}
# ```
# where $M$ is a sufficiently large negative number that is a lower bound for
# $f$ over the domain of $x$.

# As more cuts are added:
# ```math
# \lim_{K \rightarrow \infty} f^K = \min\limits_{x} f(x)
# ```

# #### Bounds

# By convexity, $f^K \le f(x)$ for all $x$. Thus, if $x^*$ is a minimizer of
# $f$, then at any point in time we can construct a lower bound for $f(x^*)$ by
# solving $f^K$.

# Moreover, since any feasible point is an upper bound, we can use the primal
# solution $x^K$ returned by solving $f^K$ to evaluate $f(x_K)$ to generate an
# upper bound.

# Therefore, $f(x^*) \in [f^K, f(x_K)]$.

# #### Implementation

# Here is pseudo-code fo the algorithm:

# 1. Take as input a function $f$ and a iteration limit $K_{max}$. Set $K = 0$,
#    and initialize $f^K$
# 2. Solve $f^K$ to obtain a candidate solution $x_{K+1}$.
# 3. Add a cut $\theta \ge f(x_{K+1}) + \frac{df}{dx}\left(x_{K+1}\right)^\top (x - x_{K+1})$ to form $f^{K+1}$.
# 4. Increment $K$
# 5. If $K = K_{max}$ STOP, otherwise, go to step 2.

# And here's a complete implementation:

function kelleys_cutting_plane(
    ## The function to be minimized.
    f::Function,
    ## The gradient of `f`. By default, we use automatic differentiation to
    ## compute the gradient of f so the user doesn't have to!
    dfdx::Function = x -> ForwardDiff.gradient(f, x);
    ## The number of arguments to `f`.
    input_dimension::Int,
    ## A lower bound for the function `f` over its domain.
    lower_bound::Float64,
    ## The number of iterations to run Kelley's algorithm for before stopping.
    iteration_limit::Int,
)
    ## Step (1):
    K = 0
    model = JuMP.Model(GLPK.Optimizer)
    JuMP.@variable(model, θ >= lower_bound)
    JuMP.@variable(model, x[1:input_dimension])
    JuMP.@objective(model, Min, θ)
    while true
        ## Step (2)
        JuMP.optimize!(model)
        x_k = JuMP.value.(x)
        ## Step (3):
        c = JuMP.@constraint(model, θ >= f(x_k) + dfdx(x_k)' * (x .- x_k))
        ## Step (4):
        K = K + 1
        ## Step (5):
        if K == iteration_limit
            break
        end
    end
    θ_K, x_K = JuMP.value(θ), JuMP.value.(x)
    println("Found solution:")
    println("  x_K   = ", x_K)
    println("  f(x*) ∈ [", θ_K, ", ", f(x_K), "]")
    return
end

# Let's run our algorithm to see what happens:

kelleys_cutting_plane(
    input_dimension = 2,
    lower_bound = 0.0,
    iteration_limit = 20,
) do x
    return (x[1] - 1)^2 + (x[2] + 2)^2
end

# ### Approximating the cost-to-go term

# In [Basic I: first steps](@ref), we discussed how you could formulate an
# optimal policy to a multistage stochastic program using the dynamic
# programming recursion:
# ```math
# \begin{aligned}
# V_i(x, \omega) = \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) + \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}[V_j(x^\prime, \varphi)]\\
# & x^\prime = T_i(\bar{x}, u, \omega) \\
# & u \in U_i(\bar{x}, \omega) \\
# & \bar{x} = x
# \end{aligned}
# ```
# where our decision rule, $\pi_i(x, \omega)$, solves this optimization problem
# and returns a $u^*$ corresponding to an optimal solution. Moreover, we alluded
# to the fact that the cost-to-go term (the nasty recursive expectation) makes
# this problem intractable to solve.

# However, if, excluding the cost-to-go term, $V_i(x, \omega)$ can be formulated
# as a linear program (this also works for convex programs, but the math is more
# involved), then we can make some progress.

# First, notice that $x$ only appears as a right-hand side term of $V_i$.
# Therefore, $V_i(x, \cdot)$ is convex with respect to $x$ for fixed $\omega$.
# Moreover, the reduced cost of the decision variable $\bar{x}$ is a subgradient
# of the function $V_i$ with respect to $x$! (This is one reason why we add the
# $\bar{x}$ and the fishing constraint $\bar{x} = x$.)

# Second, a convex combination of convex functions is also convex, so the
# cost-to-go term is a convex function of $x^\prime$.

# Stochastic dual dynamic programming converts this problem into a tractable
# form by applying Kelley's cutting plane algorithm to the cost-to-go term:
# ```math
# \begin{aligned}
# V_i^K(x, \omega) = \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) + \theta\\
# & x^\prime = T_i(\bar{x}, u, \omega) \\
# & u \in U_i(\bar{x}, \omega) \\
# & \bar{x} = x \\
# & \theta \ge \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi) + \frac{dV_j^k}{dx^\prime}\left(x^\prime_k, \varphi\right)^\top (x^\prime - x^\prime_k)\right],\quad k=1,\ldots,K \\
# & \theta \ge M
# \end{aligned}
# ```

# All we need now is a way of generating these cutting planes in an iterative
# manner.

# ### The forward pass

# ### The backward pass

# ## Implementation

# For this implementation of SDDP, we're going to try and keep things as simple
# as possible. This is very much a "vanilla" version of SDDP; it doesn't have
# (m)any fancy computational tricks that you need to code a performant or stable
# version that will work on realistic instances.

# In the interests of brevity, we will also include minimal error checking.
# Think about all the different ways you could break this code!

# ### Structs

# The first struct we are going to use is a `State` struct that will wrap an
# incoming and outgoing state variable.

struct State
    in::JuMP.VariableRef
    out::JuMP.VariableRef
end

# Next, we need a struct to wrap all of the uncertainty within a node.

struct Uncertainty
    parameterize::Function
    Ω::Vector{Any}
    P::Vector{Float64}
end

# `parameterize` is a function, which takes a realization of the random variable
# $\omega\in\Omega$ and updates the subproblem accordingly. The finite discrete
# random variable is defined by the vectors `Ω` and `P`, so that the random
# variable takes the value `Ω[i]` with probability `P[i]`. As such, `P` should
# sum to 1. (We don't check this here, but we should; we do in SDDP.jl.)

# It's also going to be useful to have a function that samples a realization of
# the random variable defined by `Ω` and `P`:

function sample_uncertainty(uncertainty::Uncertainty)
    r = rand()
    for (p, ω) in zip(uncertainty.P, uncertainty.Ω)
        if r <= p
            return ω
        end
        r -= p
    end
    error("We should never get here because P should sum to 1.0.")
end

# You should be able to work out what is going on. `rand()` samples a uniform
# random variable in `[0, 1)`.

# Now we have two building blocks, we can declare the structure of each node.

struct Node
    subproblem::JuMP.Model
    states::Dict{Symbol, State}
    uncertainty::Uncertainty
    cost_to_go::JuMP.VariableRef
end

# `subproblem` is going to be the JuMP model that we build at each node.
# `states` is a dictionary that maps a symbolic name of a state variable to a
# `State` object wrapping the incoming and outgoing state variables in
# `subproblem`. `uncertainty` is an `Uncertainty` object described above, and
# `cost_to_go` is a JuMP variable that approximates the cost-to-go term.

# Finally, out simplified policy graph is just a vector of nodes.

struct LinearPolicyGraph
    nodes::Vector{Node}
end

# We also define a nice `show` method so that we don't accidentally print a
# large amount of information to the screen.

function Base.show(io::IO, model::LinearPolicyGraph)
    return print(io, "A policy graph with $(length(model.nodes)) nodes")
end

# ### Interface functions

# Now we have some basic types, let's implment some functions so that the user
# can create a model.

# First, we need an exmaple of a function that the user will provide. Like
# SDDP.jl, this takes an empty `subproblem`, and a node index, in this case
# `t::Int`. You could change this function to change the model, or define a new
# one later in the code.

function subproblem_builder(subproblem::JuMP.Model, t::Int)
    ## Define the state variables. Note how we fix the incoming state to the
    ## initial state variable regardless of `t`! This isn't strictly necessary;
    ## it only matters that we do it for the first node.
    JuMP.@variable(subproblem, volume_in == 200)
    JuMP.@variable(subproblem, 0 <= volume_out <= 200)
    states = Dict(:volume => State(volume_in, volume_out))
    ## Define the control variables.
    JuMP.@variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        inflow
    end)
    ## Define the constraints
    JuMP.@constraints(subproblem, begin
        volume_out == volume_in + inflow - hydro_generation - hydro_spill
        demand_constraint, thermal_generation + hydro_generation == 150.0
    end)
    ## Define the objective for each stage `t`. Note that we can use `t` as an
    ## index for t = 1, 2, 3.
    fuel_cost = [50.0, 100.0, 150.0]
    JuMP.@objective(subproblem, Min, fuel_cost[t] * thermal_generation)
    ## Finally, we define the uncertainty object. Because this is a simplified
    ## implementation of SDDP, we shall politely ask the user to only modify the
    ## constraints, and not the objective function! (Not that it changes the
    ## algorithm, we just have to add more information to keep track of things.)
    uncertainty = Uncertainty([0.0, 50.0, 100.0], [1 / 3, 1 / 3, 1 / 3]) do ω
        JuMP.fix(inflow, ω)
    end
    return states, uncertainty
end

# If you've read [Basic II: adding uncertainty](@ref), this example should be
# familiar. You can probably see how some of the SDDP.jl functionality like
# [`@stageobjective`](@ref) and [`SDDP.parameterize`](@ref) help smooth some of
# the usability issues like needing to construct both the incoming and outgoing
# state variables, or needing to explicitly `return states, uncertainty`.

# The next function we need to define is the analog of
# [`SDDP.LinearPolicyGraph`](@ref). It should be pretty readable.

function LinearPolicyGraph(
    subproblem_builder::Function;
    stages::Int,
    lower_bound::Float64,
    optimizer,
)
    nodes = Node[]
    for t = 1:stages
        ## Create a model.
        model = JuMP.Model(optimizer)
        ## Use the provided function to build out each subproblem. The user's
        ## function returns a dictionary mapping `Symbol`s to `State` objects,
        ## and an `Uncertainty` object.
        states, uncertainty = subproblem_builder(model, t)
        ## Now add the cost-to-go terms:
        JuMP.@variable(model, cost_to_go >= lower_bound)
        obj = JuMP.objective_function(model)
        JuMP.@objective(model, Min, obj + cost_to_go)
        ## In the final stage, the cost-to-go is 0.0.
        if t == stages
            JuMP.fix(cost_to_go, 0.0; force = true)
        end
        push!(nodes, Node(model, states, uncertainty, cost_to_go))
    end
    return LinearPolicyGraph(nodes)
end

# ### The forward pass

# Now we're ready to code the forward pass. It takes a `::LinearPolicyGraph`,
# and returns a tuple of two things: a vector of the outgoing state variables
# visited, and a `Float64` of the cumulative stage costs that were incurred
# along the forward pass.

function forward_pass(model::LinearPolicyGraph, io::IO = stdout)
    println(io, "| Forward Pass")
    ## First, get the value of the state at the root node (e.g., x_R).
    incoming_state = Dict(
        k => JuMP.fix_value(v.in) for (k, v) in model.nodes[1].states
    )
    ## `simulation_cost` is an accumlator that is going to sum the stage-costs
    ## incurred over the forward pass.
    simulation_cost = 0.0
    ## We also need to record the outgoing state variables so we can pass them
    ## to the backward pass.
    outgoing_states = Dict{Symbol, Float64}[]
    ## Now's the meat of the forward pass: loop through each of the nodes
    for (t, node) in enumerate(model.nodes)
        println(io, "| | Visiting node $(t)")
        ## Sample the uncertainty:
        ω = sample_uncertainty(node.uncertainty)
        println(io, "| |  ω = ", ω)
        ## Before parameterizing the subproblem using the user-provided
        ## function:
        node.uncertainty.parameterize(ω)
        println(io, "| |  x = ", incoming_state)
        ## Update the incoming state variable:
        for (k, v) in incoming_state
            JuMP.fix(node.states[k].in, v; force = true)
        end
        ## Now solve the subproblem and check we found an optimal solution:
        JuMP.optimize!(node.subproblem)
        if JuMP.termination_status(node.subproblem) != JuMP.MOI.OPTIMAL
            error("Something went terribly wrong!")
        end
        ## Compute the outgoing state variables, and save them in
        ## `outgoing_states`:
        outgoing_state = Dict(k => JuMP.value(v.out) for (k, v) in node.states)
        push!(outgoing_states, outgoing_state)
        println(io, "| |  x′ = ", outgoing_state)
        ## We also need to compute the stage cost to add to our
        ## `simulation_cost` accumulator:
        stage_cost = JuMP.objective_value(node.subproblem) - JuMP.value(node.cost_to_go)
        simulation_cost += stage_cost
        println(io, "| |  C(x, u, ω) = ", stage_cost)
        ## As a final step, set the outgoing state of stage t and the incoming
        ## state of stage t + 1:
        incoming_state = outgoing_state
    end
    return outgoing_states, simulation_cost
end

# ### The backward pass

# We're now ready to code the backward pass. This is going to take a
# `::LinearPolicyGraph` object, and a vector of outgoing states from the forward
# pass. It returns a `Float64` that is a valid lower bound for the objective of
# the multistage stochastic program.

function backward_pass(
    model::LinearPolicyGraph,
    outgoing_states::Vector{Dict{Symbol, Float64}},
)
    println("| Backward pass")
    ## For the backward pass, we walk back up the nodes, from the final
    ## node to the second (we will solve the first node after this loop).
    for t = length(outgoing_states):-1:2
        println("| | Visiting node $(t)")
        ## At each step in the backward pass, we are going to solve problems in
        ## stage t, but add cuts to stage t - 1. Thus, we need:
        node_t = model.nodes[t]
        node_t1 = model.nodes[t - 1]
        ## First, fix the incoming state variables of stage t to the value of
        ## the outgoing state variables in stage t - 1.
        for (k, v) in outgoing_states[t - 1]
            JuMP.fix(node_t.states[k].in, v; force = true)
        end
        ## Then, create an empty affine expression that we will use to build
        ## up the cut expression.
        cut_expression = JuMP.AffExpr(0.0)
        ## Now for each possible realization of the uncertainty, solve the
        ## stage t subproblem, and add `p * [y + λᵀ(x - x_k)]` to the cut
        ## expression. (See the Theory section above is this isn't obvious why.)
        for (p, ω) in zip(node_t.uncertainty.P, node_t.uncertainty.Ω)
            println("| | | Solving ω = ", ω)
            node_t.uncertainty.parameterize(ω)
            JuMP.optimize!(node_t.subproblem)
            y = JuMP.objective_value(node_t.subproblem)
            println("| | |  y = ", y)
            λ = Dict(k => JuMP.reduced_cost(v.in) for (k, v) in node_t.states)
            println("| | |  λ = ", λ)
            cut_expression += p * JuMP.@expression(
                node_t1.subproblem,
                y + sum(
                    λ[k] * (x.out - outgoing_states[t - 1][k])
                    for (k, x) in node_t1.states
                ),
            )
        end
        ## And then refine the cost-to-go variable by adding a cut that is the
        ## expectation of the cuts computed in the step above.
        c = JuMP.@constraint(
            node_t1.subproblem, node_t1.cost_to_go >= cut_expression
        )
        println("| | | Adding cut : ", c)
    end
    ## Finally, compute a lower bound for the problem by evaluating the
    ## first-stage subproblem.
    first_node = model.nodes[1]
    lower_bound = 0.0
    for (p, ω) in zip(first_node.uncertainty.P, first_node.uncertainty.Ω)
        first_node.uncertainty.parameterize(ω)
        JuMP.optimize!(first_node.subproblem)
        lower_bound += p * JuMP.objective_value(first_node.subproblem)
    end
    return lower_bound
end

# Thirdly, we need a function to simulate the policy. This is going be very
# simple. It doesn't have an bells and whistles like being able to record the
# control variables.

function simulate(model::LinearPolicyGraph; replications::Int)
    ## Pipe the output to `devnull` so we don't print too much!
    simulations = [forward_pass(model, devnull) for _ = 1:replications]
    z = [s[2] for s in simulations]
    μ  = Statistics.mean(z)
    tσ = 1.96 * Statistics.std(z) / sqrt(replications)
    println("Upper bound = $(μ) ± $(tσ)")
    return simulations
end

# Finally, the `train` loop of SDDP just applies the forward and backward passes
# iteratively, followed by a final simulation to compute the upper bound
# confidence interval.

function train(
    model::LinearPolicyGraph;
    iteration_limit::Int,
    replications::Int,
)
    for i = 1:iteration_limit
        println("Starting iteration $(i)")
        outgoing_states, simulation = forward_pass(model)
        lower_bound = backward_pass(model, outgoing_states)
        println("| Finished iteration")
        println("| | simulation = ", simulation)
        println("| | lower_bound = ", lower_bound)
    end
    simulate(model; replications = replications)
    return
end

# ### Example

# First, create the model using the `subproblem_builder` function we defined
# earlier:

model = LinearPolicyGraph(
    subproblem_builder;
    stages = 3,
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
)

# Then, train a policy:

train(model; iteration_limit = 3, replications = 100)

# Success! We solved a multistage stochastic program using stochastic dual
# dynamic programming.
