# # Introductory theory

# !!! note
#     This tutorial is aimed at advanced undergraduates or early-stage graduate
#     students. You don't need prior exposure to stochastic programming!
#     (Indeed, it may be better if you don't, because our approach is
#     non-standard in the literature.)
#
#     This tutorial is also a living document. If parts are unclear, please
#     [open an issue](https://github.com/odow/SDDP.jl/issues/new) so it can be
#     improved!

# This tutorial will teach you how the stochastic dual dynamic programming
# algorithm works by implementing a simplified version of the algorithm.

# Our implementation is very much a "vanilla" version of SDDP; it doesn't have
# (m)any fancy computational tricks (e.g., the ones included in SDDP.jl) that
# you need to code a performant or stable version that will work on realistic
# instances. However, our simplified implementation will work on arbitrary
# policy graphs, including those with cycles such as infinite horizon problems!

# **Packages**
#
# This tutorial uses the following packages. For clarity, we call
# `import PackageName` so that we must prefix `PackageName.` to all functions
# and structs provided by that package. Everything not prefixed is either part
# of base Julia, or we wrote it.

import ForwardDiff
import GLPK
import JuMP
import Statistics

# !!! tip
#     You can follow along by installing the above packages, and copy-pasting
#     the code we will write into a Julia REPL. Alternatively, you can download
#     the Julia `.jl` file which created this tutorial [from Github](https://github.com/odow/SDDP.jl/blob/master/docs/src/tutorial/21_theory_intro.jl).

# ## Preliminaries: background theory

# Start this tutorial by reading [An introduction to SDDP.jl](@ref), which
# introduces the necessary notation and vocabulary that we need for this
# tutorial.

# ## Preliminaries: Kelley's cutting plane algorithm

# Kelley's cutting plane algorithm is an iterative method for minimizing convex
# functions. Given a convex function $f(x)$, Kelley's constructs an
# under-approximation of the function at the minimum by a set of first-order
# Taylor series approximations (called **cuts**) constructed at a set of points
# $k = 1,\ldots,K$:
# ```math
# \begin{aligned}
# f^K = \min\limits_{\theta \in \mathbb{R}, x \in \mathbb{R}^N} \;\; & \theta\\
# & \theta \ge f(x_k) + \frac{d}{dx}f(x_k)^\top (x - x_k),\quad k=1,\ldots,K\\
# & \theta \ge M,
# \end{aligned}
# ```
# where $M$ is a sufficiently large negative number that is a lower bound for
# $f$ over the domain of $x$.

# Kelley's cutting plane algorithm is a structured way of choosing points $x_k$
# to visit, so that as more cuts are added:
# ```math
# \lim_{K \rightarrow \infty} f^K = \min\limits_{x \in \mathbb{R}^N} f(x)
# ```
# However, before we introduce the algorithm, we need to introduce some bounds.

# ### Bounds

# By convexity, $f^K \le f(x)$ for all $x$. Thus, if $x^*$ is a minimizer of
# $f$, then at any point in time we can construct a lower bound for $f(x^*)$ by
# solving $f^K$.

# Moreover, we can use the primal solutions $x_k^*$ returned by solving $f^k$ to
# evaluate $f(x_k^*)$ to generate an upper bound.

# Therefore, $f^K \le f(x^*) \le \min\limits_{k=1,\ldots,K} f(x_k^*)$.

# When the lower bound is sufficiently close to the upper bound, we can
# terminate the algorithm and declare that we have found an solution that is
# close to optimal.

# ### Implementation

# Here is pseudo-code fo the Kelley algorithm:

# 1. Take as input a convex function $f(x)$ and a iteration limit $K_{max}$.
#    Set $K = 0$, and initialize $f^K$. Set $lb = -\infty$ and $ub = \infty$.
# 2. Solve $f^K$ to obtain a candidate solution $x_{K+1}$.
# 3. Update $lb = f^K$ and $ub = \min\{ub, f(x_{K+1})\}$.
# 4. Add a cut $\theta \ge f(x_{K+1}) + \frac{d}{dx}f\left(x_{K+1}\right)^\top (x - x_{K+1})$ to form $f^{K+1}$.
# 5. Increment $K$.
# 6. If $K = K_{max}$ or $|ub - lb| < \epsilon$, STOP, otherwise, go to step 2.

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
    ## The absolute tolerance ϵ to use for convergence.
    tolerance::Float64 = 1e-6,
)
    ## Step (1):
    K = 0
    model = JuMP.Model(GLPK.Optimizer)
    JuMP.@variable(model, θ >= lower_bound)
    JuMP.@variable(model, x[1:input_dimension])
    JuMP.@objective(model, Min, θ)
    x_k = fill(NaN, input_dimension)
    lower_bound, upper_bound = -Inf, Inf
    while true
        ## Step (2):
        JuMP.optimize!(model)
        x_k .= JuMP.value.(x)
        ## Step (3):
        lower_bound = JuMP.objective_value(model)
        upper_bound = min(upper_bound, f(x_k))
        println("K = $K : $(lower_bound) <= f(x*) <= $(upper_bound)")
        ## Step (4):
        JuMP.@constraint(model, θ >= f(x_k) + dfdx(x_k)' * (x .- x_k))
        ## Step (5):
        K = K + 1
        ## Step (6):
        if K == iteration_limit
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

kelleys_cutting_plane(
    input_dimension = 2,
    lower_bound = 0.0,
    iteration_limit = 20,
) do x
    return (x[1] - 1)^2 + (x[2] + 2)^2 + 1.0
end

# !!! warning
#     It's hard to choose a valid lower bound! If you choose one too loose, the
#     algorithm can take a long time to converge. However, if you choose one so
#     tight that $M > f(x^*)$, then you can obtain a suboptimal solution. For a
#     deeper discussion of the implications for SDDP.jl, see
#     [Choosing an initial bound](@ref).

# ## Preliminaries: approximating the cost-to-go term

# In the background theory section, we discussed how you could formulate an
# optimal policy to a multistage stochastic program using the dynamic
# programming recursion:
# ```math
# \begin{aligned}
# V_i(x, \omega) = \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) + \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}[V_j(x^\prime, \varphi)]\\
# & x^\prime = T_i(\bar{x}, u, \omega) \\
# & u \in U_i(\bar{x}, \omega) \\
# & \bar{x} = x,
# \end{aligned}
# ```
# where our decision rule, $\pi_i(x, \omega)$, solves this optimization problem
# and returns a $u^*$ corresponding to an optimal solution. Moreover, we alluded
# to the fact that the cost-to-go term (the nasty recursive expectation) makes
# this problem intractable to solve.

# However, if, excluding the cost-to-go term (i.e., the `SP` formulation),
# $V_i(x, \omega)$ can be formulated as a linear program (this also works for
# convex programs, but the math is more involved), then we can make some
# progress by noticing that $x$ only appears as a right-hand side term of the
# fishing constraint $\bar{x} = x$.

# Therefore, $V_i(x, \cdot)$ is convex with respect to $x$ for fixed $\omega$.
# Moreover, if we implement the constraint $\bar{x} = x$ by setting the lower-
# and upper bounds of $\bar{x}$ to $x$, then the reduced cost of the decision
# variable $\bar{x}$ is a subgradient of the function $V_i$ with respect to $x$!
# (This is the algorithmic simplification that leads us to add $\bar{x}$ and the
# fishing constraint $\bar{x} = x$.)

# !!! tip
#     The subproblem can have binary and integer variables, but you'll need to
#     use Lagrangian duality to compute a subgradient!

# Stochastic dual dynamic programming converts this problem into a tractable
# form by applying Kelley's cutting plane algorithm to the $V_j$ functions in
# the cost-to-go term:
# ```math
# \begin{aligned}
# V_i^K(x, \omega) = \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) + \theta\\
# & x^\prime = T_i(\bar{x}, u, \omega) \\
# & u \in U_i(\bar{x}, \omega) \\
# & \bar{x} = x \\
# & \theta \ge \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi) + \frac{d}{dx^\prime}V_j^k(x^\prime_k, \varphi)^\top (x^\prime - x^\prime_k)\right],\quad k=1,\ldots,K \\
# & \theta \ge M.
# \end{aligned}
# ```

# All we need now is a way of generating these cutting planes in an iterative
# manner. Before we get to that though, let's start writing some code.

# ## Implementation: modeling

# Let's make a start by defining the problem structure. Like SDDP.jl, we need a
# few things:
#
# 1. A description of the structure of the policy graph: how many nodes there
#    are, and the arcs linking the nodes together with their corresponding
#    probabilities.
# 2. A JuMP model for each node in the policy graph.
# 3. A way to identify the incoming and outgoing state variables of each node.
# 4. A description of the random variable, as well as a function that we can
#    call that will modify the JuMP model to reflect the realization of the
#    random variable.
# 5. A decision variable to act as the approximated cost-to-go term.

# !!! warning
#     In the interests of brevity, there is minimal error checking. Think about
#     all the different ways you could break the code!

# ### Structs

# The first struct we are going to use is a `State` struct that will wrap an
# incoming and outgoing state variable:

struct State
    in::JuMP.VariableRef
    out::JuMP.VariableRef
end

# Next, we need a struct to wrap all of the uncertainty within a node:

struct Uncertainty
    parameterize::Function
    Ω::Vector{Any}
    P::Vector{Float64}
end

# `parameterize` is a function which takes a realization of the random variable
# $\omega\in\Omega$ and updates the subproblem accordingly. The finite discrete
# random variable is defined by the vectors `Ω` and `P`, so that the random
# variable takes the value `Ω[i]` with probability `P[i]`. As such, `P` should
# sum to 1. (We don't check this here, but we should; we do in SDDP.jl.)

# Now we have two building blocks, we can declare the structure of each node:

struct Node
    subproblem::JuMP.Model
    states::Dict{Symbol,State}
    uncertainty::Uncertainty
    cost_to_go::JuMP.VariableRef
end

# * `subproblem` is going to be the JuMP model that we build at each node.
# * `states` is a dictionary that maps a symbolic name of a state variable to a
#   `State` object wrapping the incoming and outgoing state variables in
#   `subproblem`.
# * `uncertainty` is an `Uncertainty` object described above.
# * `cost_to_go` is a JuMP variable that approximates the cost-to-go term.

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
# large amount of information to the screen when creating a model:

function Base.show(io::IO, model::PolicyGraph)
    println(io, "A policy graph with $(length(model.nodes)) nodes")
    println(io, "Arcs:")
    for (from, arcs) in enumerate(model.arcs)
        for (to, probability) in arcs
            println(io, "  $(from) => $(to) w.p. $(probability)")
        end
    end
    return
end

# ### Functions

# Now we have some basic types, let's implement some functions so that the user
# can create a model.

# First, we need an example of a function that the user will provide. Like
# SDDP.jl, this takes an empty `subproblem`, and a node index, in this case
# `t::Int`. You could change this function to change the model, or define a new
# one later in the code.

# We're going to copy the example from [An introduction to SDDP.jl](@ref),
# with some minor adjustments for the fact we don't have many of the bells and
# whistles of SDDP.jl. You can probably see how some of the SDDP.jl
# functionality like [`@stageobjective`](@ref) and [`SDDP.parameterize`](@ref)
# help smooth some of the usability issues like needing to construct both the
# incoming and outgoing state variables, or needing to explicitly declare
# `return states, uncertainty`.

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
        hydro_generation >= 0
        hydro_spill >= 0
        inflow
    end)
    ## Define the constraints
    JuMP.@constraints(
        subproblem,
        begin
            volume_out == volume_in + inflow - hydro_generation - hydro_spill
            demand_constraint, thermal_generation + hydro_generation == 150.0
        end
    )
    ## Define the objective for each stage `t`. Note that we can use `t` as an
    ## index for t = 1, 2, 3.
    fuel_cost = [50.0, 100.0, 150.0]
    JuMP.@objective(subproblem, Min, fuel_cost[t] * thermal_generation)
    ## Finally, we define the uncertainty object. Because this is a simplified
    ## implementation of SDDP, we shall politely ask the user to only modify the
    ## constraints, and not the objective function! (Not that it changes the
    ## algorithm, we just have to add more information to keep track of things.)
    uncertainty = Uncertainty([0.0, 50.0, 100.0], [1 / 3, 1 / 3, 1 / 3]) do ω
        return JuMP.fix(inflow, ω)
    end
    return states, uncertainty
end

# The next function we need to define is the analog of
# [`SDDP.PolicyGraph`](@ref). It should be pretty readable.

function PolicyGraph(
    subproblem_builder::Function;
    graph::Vector{Dict{Int,Float64}},
    lower_bound::Float64,
    optimizer,
)
    nodes = Node[]
    for t in 1:length(graph)
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

# Then, we can create a model using the `subproblem_builder` function we defined
# earlier:

model = PolicyGraph(
    subproblem_builder;
    graph = [Dict(2 => 1.0), Dict(3 => 1.0), Dict{Int,Float64}()],
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
)

# ## Implementation: helpful samplers

# Before we get properly coding the solution algorithm, it's also going to be
# useful to have a function that samples a realization of the random variable
# defined by `Ω` and `P`.

function sample_uncertainty(uncertainty::Uncertainty)
    r = rand()
    for (p, ω) in zip(uncertainty.P, uncertainty.Ω)
        r -= p
        if r < 0.0
            return ω
        end
    end
    return error("We should never get here because P should sum to 1.0.")
end

# !!! note
#     `rand()` samples a uniform random variable in `[0, 1)`.

# For example:

for i in 1:3
    println("ω = ", sample_uncertainty(model.nodes[1].uncertainty))
end

# It's also going to be useful to define a function that generates a random walk
# through the nodes of the graph:

function sample_next_node(model::PolicyGraph, current::Int)
    if length(model.arcs[current]) == 0
        ## No outgoing arcs!
        return nothing
    else
        r = rand()
        for (to, probability) in model.arcs[current]
            r -= probability
            if r < 0.0
                return to
            end
        end
        ## We looped through the outgoing arcs and still have probability left
        ## over! This means we've hit an implicit "zero" node.
        return nothing
    end
end

# For example:

for i in 1:3
    ## We use `repr` to print the next node, because `sample_next_node` can
    ## return `nothing`.
    println("Next node from $(i) = ", repr(sample_next_node(model, i)))
end

# This is a little boring, because our graph is simple. However, more
# complicated graphs will generate more interesting trajectories!

# ## Implementation: the forward pass

# Recall that, after approximating the cost-to-go term, we need a way of
# generating the cuts. As the first step, we need a way of generating candidate
# solutions $x_k^\prime$. However, unlike the Kelley's example, our functions
# $V_j^k(x^\prime, \varphi)$ need two inputs: an outgoing state variable and a
# realization of the random variable.

# One way of getting these inputs is just to pick a random (feasible) value.
# However, in doing so, we might pick outgoing state variables that we will
# never see in practice, or we might infrequently pick outgoing state variables
# that we will often see in practice. Therefore, a better way of generating the
# inputs is to use a simulation of the policy, which we call the **forward**
# **pass**.

# The forward pass walks the policy graph from start to end, transitioning
# randomly along the arcs. At each node, it observes a realization of the random
# variable and solves the approximated subproblem to generate a candidate
# outgoing state variable $x_k^\prime$. The outgoing state variable is passed as
# the incoming state variable to the next node in the trajectory.

function forward_pass(model::PolicyGraph, io::IO = stdout)
    println(io, "| Forward Pass")
    ## First, get the value of the state at the root node (e.g., x_R).
    incoming_state =
        Dict(k => JuMP.fix_value(v.in) for (k, v) in model.nodes[1].states)
    ## `simulation_cost` is an accumlator that is going to sum the stage-costs
    ## incurred over the forward pass.
    simulation_cost = 0.0
    ## We also need to record the nodes visited and resultant outgoing state
    ## variables so we can pass them to the backward pass.
    trajectory = Tuple{Int,Dict{Symbol,Float64}}[]
    ## Now's the meat of the forward pass: beginning at the first node:
    t = 1
    while t !== nothing
        node = model.nodes[t]
        println(io, "| | Visiting node $(t)")
        ## Sample the uncertainty:
        ω = sample_uncertainty(node.uncertainty)
        println(io, "| | | ω = ", ω)
        ## Parameterizing the subproblem using the user-provided function:
        node.uncertainty.parameterize(ω)
        println(io, "| | | x = ", incoming_state)
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
        println(io, "| | | x′ = ", outgoing_state)
        ## We also need to compute the stage cost to add to our
        ## `simulation_cost` accumulator:
        stage_cost =
            JuMP.objective_value(node.subproblem) - JuMP.value(node.cost_to_go)
        simulation_cost += stage_cost
        println(io, "| | | C(x, u, ω) = ", stage_cost)
        ## As a penultimate step, set the outgoing state of stage t and the
        ## incoming state of stage t + 1, and add the node to the trajectory.
        incoming_state = outgoing_state
        push!(trajectory, (t, outgoing_state))
        ## Finally, sample a new node to step to. If `t === nothing`, the
        ## `while` loop will break.
        t = sample_next_node(model, t)
    end
    return trajectory, simulation_cost
end

# Let's take a look at one forward pass:

trajectory, simulation_cost = forward_pass(model);

# ## Implementation: the backward pass

# From the forward pass, we obtained a vector of nodes visted and their
# corresponding outgoing state variables. Now we need to refine the
# approximation for each node at the candidate solution for the outgoing state
# variable. That is, we need to add a new cut:
# ```math
# \theta \ge \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi) + \frac{d}{dx^\prime}V_j^k(x^\prime_k, \varphi)^\top (x^\prime - x^\prime_k)\right]
# ```
# or alternatively:
# ```math
# \theta \ge \sum\limits_{j \in i^+} \sum\limits_{\varphi \in \Omega_j} p_{ij} p_{\varphi}\left[V_j^k(x^\prime_k, \varphi) + \frac{d}{dx^\prime}V_j^k(x^\prime_k, \varphi)^\top (x^\prime - x^\prime_k)\right]
# ```

# It doesn't matter what order we visit the nodes to generate these cuts for.
# For example, we could compute them all in parallel, using the current
# approximations of $V^K_i$.

# However, we can be smarter than that.

# If we traverse the list of nodes visited in the forward pass in reverse, then
# we come to refine the $i$th node in the trajectory, we will already have
# improved the approximation of the $(i+1)$th node in the trajectory as well!
# Therefore, our refinement of the $i$th node will be better than if we improved
# node $i$ first, and then refined node $(i+1)$.

# Because we walk the nodes in reverse, we call this the **backward pass**.

# !!! info
#     If you're into deep learning, you could view this as the equivalent of
#     back-propagation: the forward pass pushes primal information through the
#     graph (outgoing state variables), and the backward pass pulls dual
#     information (cuts) back through the graph to improve our decisions on the
#     next forward pass.

function backward_pass(
    model::PolicyGraph,
    trajectory::Vector{Tuple{Int,Dict{Symbol,Float64}}},
    io::IO = stdout,
)
    println(io, "| Backward pass")
    ## For the backward pass, we walk back up the nodes.
    for i in reverse(1:length(trajectory))
        index, outgoing_states = trajectory[i]
        node = model.nodes[index]
        println(io, "| | Visiting node $(index)")
        if length(model.arcs[index]) == 0
            ## If there are no children, the cost-to-go is 0.
            println(io, "| | | Skipping node because the cost-to-go is 0")
            continue
        end
        ## Create an empty affine expression that we will use to build up the
        ## right-hand side of the cut expression.
        cut_expression = JuMP.AffExpr(0.0)
        ## For each node j ∈ i⁺
        for (j, P_ij) in model.arcs[index]
            next_node = model.nodes[j]
            ## Set the incoming state variables of node j to the outgoing state
            ## variables of node i
            for (k, v) in outgoing_states
                JuMP.fix(next_node.states[k].in, v; force = true)
            end
            ## Then for each realization of φ ∈ Ωⱼ
            for (pφ, φ) in zip(next_node.uncertainty.P, next_node.uncertainty.Ω)
                ## Setup and solve for the realization of φ
                println(io, "| | | Solving φ = ", φ)
                next_node.uncertainty.parameterize(φ)
                JuMP.optimize!(next_node.subproblem)
                ## Then prepare the cut `P_ij * pφ * [V + dVdxᵀ(x - x_k)]``
                V = JuMP.objective_value(next_node.subproblem)
                println(io, "| | | | V = ", V)
                dVdx = Dict(
                    k => JuMP.reduced_cost(v.in) for (k, v) in next_node.states
                )
                println(io, "| | | | dVdx′ = ", dVdx)
                cut_expression += JuMP.@expression(
                    node.subproblem,
                    P_ij *
                    pφ *
                    (
                        V + sum(
                            dVdx[k] * (x.out - outgoing_states[k]) for
                            (k, x) in node.states
                        )
                    ),
                )
            end
        end
        ## And then refine the cost-to-go variable by adding the cut:
        c = JuMP.@constraint(node.subproblem, node.cost_to_go >= cut_expression)
        println(io, "| | | Adding cut : ", c)
    end
    return nothing
end

# ## Implementation: bounds

# ### Lower bounds

# Recall from Kelley's that we can obtain a lower bound for $f(x^*)$ be
# evaluating $f^K$. The analogous lower bound for a multistage stochastic
# program is:

# ```math
# \mathbb{E}_{i \in R^+, \omega \in \Omega_i}[V_i^K(x_R, \omega)] \le \min_{\pi} \mathbb{E}_{i \in R^+, \omega \in \Omega_i}[V_i^\pi(x_R, \omega)]
# ```

# Here's how we compute the lower bound:

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

# !!! note
#     The implementation is simplified because we assumed that there is only one
#     arc from the root node, and that it pointed to the first node in the
#     vector.

# Because we haven't trained a policy yet, the lower bound is going to be very
# bad:

lower_bound(model)

# ### Upper bounds

# With Kelley's algorithm, we could easily construct an upper bound by
# evaluating $f(x_K)$. However, it is almost always intractable to evaluate an
# upper bound for multistage stochastic programs due to the large number of
# nodes and the nested expectations. Instead, we can perform a Monte Carlo
# simulation of the policy to build a statistical estimate for the value of
# $\mathbb{E}_{i \in R^+, \omega \in \Omega_i}[V_i^\pi(x_R, \omega)]$, where
# $\pi$ is the policy defined by the current approximations $V^K_i$.

function upper_bound(model::PolicyGraph; replications::Int)
    ## Pipe the output to `devnull` so we don't print too much!
    simulations = [forward_pass(model, devnull) for i in 1:replications]
    z = [s[2] for s in simulations]
    μ = Statistics.mean(z)
    tσ = 1.96 * Statistics.std(z) / sqrt(replications)
    return μ, tσ
end

# !!! note
#     The width of the confidence interval is incorrect if there are cycles in
#     the graph, because the distribution of simulation costs `z` is not
#     symmetric. The mean is correct, however.

# ### Termination criteria

# In Kelley's algorithm, the upper bound was deterministic. Therefore, we could
# terminate the algorithm when the lower bound was sufficiently close to the
# upper bound. However, our upper bound for SDDP is not deterministic; it is a
# confidence interval!

# Some people suggest terminating SDDP when the lower bound is contained within
# the confidence interval. However, this is a poor choice because it is too easy
# to generate a false positive. For example, if we use a small number of
# replications then the width of the confidence will be large, and we are more
# likely to terminate!

# In a future tutorial (not yet written...) we will discuss termination criteria
# in more depth. For now, pick a large number of iterations and train for as
# long as possible.

# !!! tip
#     For a rule of thumb, pick a large number of iterations to train the
#     policy for (e.g.,
#     $10 \times |\mathcal{N}| \times \max\limits_{i\in\mathcal{N}} |\Omega_i|$)

# ## Implementation: the training loop

# The `train` loop of SDDP just applies the forward and backward passes
# iteratively, followed by a final simulation to compute the upper bound
# confidence interval:

function train(
    model::PolicyGraph;
    iteration_limit::Int,
    replications::Int,
    io::IO = stdout,
)
    for i in 1:iteration_limit
        println(io, "Starting iteration $(i)")
        outgoing_states, _ = forward_pass(model, io)
        backward_pass(model, outgoing_states, io)
        println(io, "| Finished iteration")
        println(io, "| | lower_bound = ", lower_bound(model))
    end
    println(io, "Termination status: iteration limit")
    μ, tσ = upper_bound(model; replications = replications)
    println(io, "Upper bound = $(μ) ± $(tσ)")
    return
end

# Using our `model` we defined earlier, we can go:

train(model; iteration_limit = 3, replications = 100)

# Success! We trained a policy for a finite horizon multistage stochastic
# program using stochastic dual dynamic programming.

# ## Implementation: evaluating the policy

# A final step is the ability to evaluate the policy at a given point.

function evaluate_policy(
    model::PolicyGraph;
    node::Int,
    incoming_state::Dict{Symbol,Float64},
    random_variable,
)
    the_node = model.nodes[node]
    the_node.uncertainty.parameterize(random_variable)
    for (k, v) in incoming_state
        JuMP.fix(the_node.states[k].in, v; force = true)
    end
    JuMP.optimize!(the_node.subproblem)
    return Dict(
        k => JuMP.value.(v) for
        (k, v) in JuMP.object_dictionary(the_node.subproblem)
    )
end

evaluate_policy(
    model;
    node = 1,
    incoming_state = Dict(:volume => 150.0),
    random_variable = 75,
)

# !!! note
#     The random variable can be **out-of-sample**, i.e., it doesn't have to be
#     in the vector $\Omega$ we created when defining the model! This is a
#     notable difference to other multistage stochastic solution methods like
#     progressive hedging or using the deterministic equivalent.

# ## Example: infinite horizon

# As promised earlier, our implementation is actually pretty general. It can
# solve any multistage stochastic (linear) program defined by a policy graph,
# including infinite horizon problems!

# Here's an example, where we have extended our earlier problem with an arc from
# node 3 to node 2 with probability 0.5. You can interpret the 0.5 as a discount
# factor.

model = PolicyGraph(
    subproblem_builder;
    graph = [Dict(2 => 1.0), Dict(3 => 1.0), Dict(2 => 0.5)],
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
)

# Then, train a policy:

train(model; iteration_limit = 3, replications = 100)

# Success! We trained a policy for an infinite horizon multistage stochastic
# program using stochastic dual dynamic programming. Note how some of the
# forward passes are different lengths!

evaluate_policy(
    model;
    node = 3,
    incoming_state = Dict(:volume => 100.0),
    random_variable = 10.0,
)
