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
# manner. Stochastic dual dynamic programming does this in two phases: a
# **forward pass** and a **backward pass**.

# ### The forward pass

# The forward pass walks the policy graph from start to end, and solves each
# approximated subproblem to generate a candidate outgoing state variable
# $x_k^\prime$ at which to generate a cut.

# ### The backward pass

# ## Implementation

# For this implementation of SDDP, we're going to try and keep things as simple.
# This is very much a "vanilla" version of SDDP; it doesn't have (m)any fancy
# computational tricks that you need to code a performant or stable version that
# will work on realistic instances. However, it will work on arbitrary policy
# graphs, including those with cycles such as infinite horizon problems!

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

# Finally, we define a simplified policy graph as follows:
struct PolicyGraph
    nodes::Vector{Node}
    arcs::Vector{Dict{Int,Float64}}
end

# There is a vector of nodes, as well as a data structure for the arcs. `arcs`
# is a vector of dictionaries, where `arcs[i][j]` gives the probabiltiy of
# transitioning from node `i` to node `j`, if an arc exists.

# To simplify things, we will assume that the root node transitions to node `1`
# with probability 1, and there are no other incoming arcs to node 1. Notably,
# we can still define cyclic graphs though!

# We also define a nice `show` method so that we don't accidentally print a
# large amount of information to the screen.

function Base.show(io::IO, model::PolicyGraph)
    return print(io, "A policy graph with $(length(model.nodes)) nodes")
end

# It's also going to be useful to define a function that generates a random walk
# through the nodes of the graph:

function sample_graph(model::PolicyGraph)
    trajectory, current_node = Int[], 1
    finished = false
    while !finished
        push!(trajectory, current_node)
        r = rand()
        for (to, probability) in model.arcs[current_node]
            r -= probability
            if r < 0.0
                current_node = to
                break
            end
        end
        if r >= 0
            ## We looped through the outgoing arcs and still have probability
            ## left over! This means it's time to stop walking.
            finished = true
        end
    end
    return trajectory
end

# As well as a function that computes a lower bound for the objective of the
# policy graph:

function lower_bound(model::PolicyGraph)
    node = model.nodes[1]
    bound = 0.0
    for (p, ω) in zip(node.uncertainty.P, node.uncertainty.Ω)
        node.uncertainty.parameterize(ω)
        JuMP.optimize!(node.subproblem)
        bound += p * JuMP.objective_value(node.subproblem)
    end
    return bound
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
# [`SDDP.PolicyGraph`](@ref). It should be pretty readable.

function PolicyGraph(
    subproblem_builder::Function;
    graph::Vector{Dict{Int,Float64}},
    lower_bound::Float64,
    optimizer,
)
    nodes = Node[]
    for t = 1:length(graph)
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
        ## If there are no outgoing arcs, the cost-to-go is 0.0.
        if length(graph[t]) == 0
            JuMP.fix(cost_to_go, 0.0; force = true)
        end
        push!(nodes, Node(model, states, uncertainty, cost_to_go))
    end
    return PolicyGraph(nodes, graph)
end

# ### The forward pass

# Now we're ready to code the forward pass. It takes a `::PolicyGraph`,
# and returns a tuple of two things: a vector of the outgoing state variables
# visited, and a `Float64` of the cumulative stage costs that were incurred
# along the forward pass.

function forward_pass(model::PolicyGraph, io::IO = stdout)
    println(io, "| Forward Pass")
    ## First, get the value of the state at the root node (e.g., x_R).
    incoming_state = Dict(
        k => JuMP.fix_value(v.in) for (k, v) in model.nodes[1].states
    )
    ## `simulation_cost` is an accumlator that is going to sum the stage-costs
    ## incurred over the forward pass.
    simulation_cost = 0.0
    ## We also need to record the nodes visited and resultant outgoing state
    ## variables so we can pass them to the backward pass.
    trajectory = Tuple{Int,Dict{Symbol,Float64}}[]
    ## Now's the meat of the forward pass: loop through each of the nodes
    for t in sample_graph(model)
        node = model.nodes[t]
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
        ## Compute the outgoing state variables:
        outgoing_state = Dict(k => JuMP.value(v.out) for (k, v) in node.states)
        println(io, "| |  x′ = ", outgoing_state)
        ## We also need to compute the stage cost to add to our
        ## `simulation_cost` accumulator:
        stage_cost = JuMP.objective_value(node.subproblem) - JuMP.value(node.cost_to_go)
        simulation_cost += stage_cost
        println(io, "| |  C(x, u, ω) = ", stage_cost)
        ## As a final step, set the outgoing state of stage t and the incoming
        ## state of stage t + 1, and add the node to the trajectory.
        incoming_state = outgoing_state
        push!(trajectory, (t, outgoing_state))
    end
    return trajectory, simulation_cost
end

# ### The backward pass

# We're now ready to code the backward pass. This is going to take a
# `::PolicyGraph` object, and a vector of (node, outgoing states) tuples from
# the forward pass. It returns a `Float64` that is a valid lower bound for the
# objective of the multistage stochastic program.

function backward_pass(
    model::PolicyGraph,
    trajectory::Vector{Tuple{Int,Dict{Symbol,Float64}}},
    io::IO = stdout,
)
    println(io, "| Backward pass")
    ## For the backward pass, we walk back up the nodes, from the final
    ## node to the second (we will solve the first node after this loop).
    for i = reverse(1:length(trajectory)-1)
        index, outgoing_states = trajectory[i]
        node = model.nodes[index]
        println(io, "| | Visiting node $(index)")
        ## Create an empty affine expression that we will use to build up the
        ## cut expression.
        cut_expression = JuMP.AffExpr(0.0)
        ## Now for each possible node and realization of the uncertainty, solve
        ## the subproblem, and add `P_ij * p_ω * [y + λᵀ(x - x_k)]` to the cut
        ## expression. (See the Theory section above is this isn't obvious why.)
        for (j, P_ij) in model.arcs[index]
            next_node = model.nodes[j]
            for (k, v) in outgoing_states
                JuMP.fix(next_node.states[k].in, v; force = true)
            end
            for (pω, ω) in zip(next_node.uncertainty.P, next_node.uncertainty.Ω)
                println(io, "| | | Solving ω = ", ω)
                next_node.uncertainty.parameterize(ω)
                JuMP.optimize!(next_node.subproblem)
                y = JuMP.objective_value(next_node.subproblem)
                println(io, "| | |  y = ", y)
                λ = Dict(k => JuMP.reduced_cost(v.in) for (k, v) in next_node.states)
                println(io, "| | |  λ = ", λ)
                cut_expression += P_ij * pω * JuMP.@expression(
                    node.subproblem,
                    y + sum(
                        λ[k] * (x.out - outgoing_states[k])
                        for (k, x) in node.states
                    ),
                )
            end
        end
        ## And then refine the cost-to-go variable by adding a cut that is the
        ## expectation of the cuts computed in the step above.
        c = JuMP.@constraint(
            node.subproblem, node.cost_to_go >= cut_expression
        )
        println(io, "| | | Adding cut : ", c)
    end
    ## Finally, compute a lower bound for the problem by evaluating the
    ## first-stage subproblem.
    return lower_bound(model)
end

# Thirdly, we need a function to simulate the policy. This is going be very
# simple. It doesn't have an bells and whistles like being able to record the
# control variables. The confidence interval is also incorrect if there are
# cycles in the graph, because the distribution of simulation costs `z` is not
# symmetric.

function simulate(model::PolicyGraph, io::IO = stdout; replications::Int)
    ## Pipe the output to `devnull` so we don't print too much!
    simulations = [forward_pass(model, devnull) for _ = 1:replications]
    z = [s[2] for s in simulations]
    μ  = Statistics.mean(z)
    tσ = 1.96 * Statistics.std(z) / sqrt(replications)
    println(io, "Upper bound = $(μ) ± $(tσ)")
    return simulations
end

# Finally, the `train` loop of SDDP just applies the forward and backward passes
# iteratively, followed by a final simulation to compute the upper bound
# confidence interval.

function train(
    model::PolicyGraph;
    iteration_limit::Int,
    replications::Int,
    io::IO = stdout,
)
    for i = 1:iteration_limit
        println(io, "Starting iteration $(i)")
        outgoing_states, simulation = forward_pass(model, io)
        lower_bound = backward_pass(model, outgoing_states, io)
        println(io, "| Finished iteration")
        println(io, "| | simulation = ", simulation)
        println(io, "| | lower_bound = ", lower_bound)
    end
    simulate(model, io; replications = replications)
    return
end

# ### Example: finite horizon

# First, create the model using the `subproblem_builder` function we defined
# earlier:

model = PolicyGraph(
    subproblem_builder;
    graph = [
        Dict(2 => 1.0),
        Dict(3 => 1.0),
        Dict{Int,Float64}(),
    ],
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
)

# Then, train a policy:

train(model; iteration_limit = 3, replications = 100)

# Success! We trained a policy for a finite horizon multistage stochastic
# program using stochastic dual dynamic programming.

# ### Example: infinite horizon

# First, create the model using the `subproblem_builder` function we defined
# earlier:

model = PolicyGraph(
    subproblem_builder;
    graph = [
        Dict(2 => 1.0),
        Dict(3 => 1.0),
        Dict(2 => 0.5),
    ],
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
)

# Then, train a policy:

train(model; iteration_limit = 3, replications = 100)

# Success! We trained a policy for an infinite horizon multistage stochastic
# program using stochastic dual dynamic programming. Note how some of the
# forward passes are different lengths!
