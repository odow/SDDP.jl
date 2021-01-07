# # Expert I: pedagogical SDDP

# In this tutorial we walk through a simplified implementation of stochastic
# dual dynamic programming to explain the key concepts.

# For this implementation of SDDP, we're going to try and keep things simple.
# This is very much a "vanilla" version of SDDP; it doesn't have (m)any fancy
# computational tricks that you need to code a performant or stable version that
# will work on realistic instances. However, it will work on arbitrary policy
# graphs, including those with cycles such as infinite horizon problems!

# !!! warning
#     In the interests of brevity, there is minimal error checking. Think about
#     all the different ways you could break the code!

# This tutorial uses the following packages. For clarity, we call
# `import PackageName` so that we must prefix `PackageName.` to all functions
# and structs provided by that package. Everything not prefixed is either part
# of base Julia, or we wrote it.

import ForwardDiff
import GLPK
import JuMP
import Statistics

# ## Preliminaries: background theory

# !!! tip
#     This section is copied verbatim from [Basic I: first steps](@ref). If it's
#     familiar, skip to [Preliminaries: Kelley's cutting plane algorithm](@ref).

# Multistage stochastic programming is complicated, and the literature has not
# settled upon standard naming conventions, so we must begin with some
# unavoidable theory and notation.

# A multistage stochastic program can be modeled by a **policy graph**. A policy
# graph is a graph with nodes and arcs. The simplest type of policy graph is a
# linear graph. Here's a linear graph with three nodes:

# ![Linear policy graph](../assets/stochastic_linear_policy_graph.png)

# In addition to nodes 1, 2, and 3, there is also a root node (0), and three
# arcs. Each arc has an origin node and a destination node, like `0 => 1`, and a
# corresponding probability of transitioning from the origin to the destination.
# For now, we can forget about the arc probabilities, because they are all 1.0.
# The squiggly lines coming into each node represent random variables; we'll
# talk about them below.

# We denote the set of nodes by $\mathcal{N}$, the root node by $R$, and the
# probability of transitioning from node $i$ to node $j$ by $p_{ij}$. (If no arc
# exists, then $p_{ij} = 0$). We define the set of successors of node $i$ as
# $i^+ = \{j \in N | P(i => j) > 0\}$.

# Each square node in the graph corresponds to a place at which the agent makes
# a decision, and we call moments in time at which the agent makes a decision
# **stages**. By convention, we try to draw policy graphs from left-to-right,
# with the stages as columns. There can be more than one node in a stage! Here's
# an example, taken from the paper [Dowson (2020)](https://doi.org/10.1002/net.21932):

# ![Markovian policy graph](../assets/powder_policy_graph.png)

# The columns represent time, and the rows represent different states of the
# world. In this case, the rows represent different prices that milk can be sold
# for at the end of each year. You can think of the nodes as forming a Markov
# chain, therefore, we call problems with a structure like this **Markovian**
# **policy graphs**. Moreover, note that policy graphs can have cycles! This
# allows them to model infinite horizon problems.

# A common feature of multistage stochastic optimization problems is that they
# model an agent controlling a system over time. This system can be described by
# three types of variables.

# 1. **State** variables track a property of the system over time.

#    Each node has an associated _incoming_ state variable (the value of the
#    state at the start of the node), and an _outgoing_ state variable (the
#    value of the state at the end of the node).
#
#    Examples of state variables include the volume of water in a reservoir, the
#    number of units of inventory in a warehouse, or the spatial position of a
#    moving vehicle.
#
#    Because state variables track the system over time, each node must have the
#    same set of state variables.
#
#    We denote state variables by the letter $x$ for the incoming state variable
#    and $x^\prime$ for the outgoing state variable.
#
# 2. **Control** variables are actions taken (implicitly or explicitly) by the
#    agent within a node which modify the state variables.
#
#    Examples of control variables include releases of water from the reservoir,
#    sales or purchasing decisions, and acceleration or braking of the vehicle.
#
#    Control variables are local to a node $i$, and they can differ between
#    nodes. For example, some control variables may be available within certain
#    nodes.
#
#    We denote control variables by the letter $u$.
#
# 3. **Random** variables are finite, discrete, exogenous random variables that
#    the agent observes at the start of a node, before the control variables are
#    decided.
#
#    Examples of random variables include rainfall inflow into a reservoir,
#    probalistic perishing of inventory, and steering errors in a vehicle.
#
#    Random variables are local to a node $i$, and they can differ between
#    nodes. For example, some nodes may have random variables, and some nodes
#    may not.
#
#    We denote random variables by the Greek letter $\omega$ and the sample
#    space from which they are drawn by $\Omega_i$. The probability of sampling
#    $\omega$ is denoted $p_{\omega}$ for simplicity.

# In a node $i$, the three variables are related by a **transition function**,
# which maps the incoming state, the controls, and the random variables to the
# outgoing state as follows: $x^\prime = T_i(x, u, \omega)$.

# As a result of entering a node $i$ with the incoming state $x$, observing
# random variable $\omega$, and choosing control $u$, the agent incurs a cost
# $C_i(x, u, \omega)$. (If the agent is a maximizer, this can be a profit, or a
# negative cost.) We call $C_i$ the **stage objective**.

# To choose their control variables in node $i$, the agent uses a **decision**
# **rule** $u = \pi_i(x, \omega)$, which is a function that maps the incoming
# state variable and observation of the random variable to a control $u$. This
# control must satisfy some feasibilty requirements $u \in U_i(x, \omega)$.

# The set of decision rules, with one element for each node in the policy graph,
# is called a **policy**.

# The goal of the agent is to find a policy that minimizes the expected cost of
# starting at the root node with some initial condition $x_R$, and proceeding
# from node to node along the probabilistic arcs until they reach a node with no
# outgoing arcs.

# ```math
# \min_{\pi} \mathbb{E}_{i \in R^+, \omega \in \Omega_i}[V_i^\pi(x_R, \omega)]
# ```
# where
# ```math
# V_i^\pi(x, \omega) = C_i(x, u, \omega) + \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}[V_j(x^\prime, \varphi)]
# ```
# where $u = \pi_i(x, \omega) \in U_i(x, \omega)$, and
# $x^\prime = T_i(x, u, \omega)$.

# The expectations are a bit complicated, but they are equivalent to:
# ```math
# \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}[V_j(x^\prime, \varphi)] = \sum\limits_{j \in i^+} p_{ij} \left[\sum\limits_{\varphi \in \Omega_j} p_{\varphi}\left[V_j(x^\prime, \varphi)\right]\right]
# ```

# An optimal policy is the set of decision rules that the agent can use to make
# these decisions and achieve the smallest expected cost.

# Often, computing the cost of a policy is intractable due to the large number
# of nodes or possible realizations of the random variables. Instead, we can
# evaluate the policy using a Monte Carlo simulation. Each replicate of the
# simulation starts at the root node and probabilistically walks along the arcs
# of the policy graph until it reaches a node with not outgoing arcs. The cost
# of a replicate is the sum of the costs incurred at each node that was visited.

# ### Dynamic programming and subproblems

# Now that we have formulated our problem, we need some ways of computing
# optimal decision rules. One way is to just use a heuristic like "choose a
# control randomly from the set of feasible controls." However, such a policy is
# unlikely to be optimal.

# One way of obtaining an optimal policy is to use Bellman's principle of
# optimality, a.k.a Dynamic Programming, and define a recursive **subproblem**
# as follows:
# ```math
# \begin{aligned}
# V_i(x, \omega) = \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) + \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}[V_j(x^\prime, \varphi)]\\
# & x^\prime = T_i(\bar{x}, u, \omega) \\
# & u \in U_i(\bar{x}, \omega) \\
# & \bar{x} = x
# \end{aligned}
# ```
# Our decision rule, $\pi_i(x, \omega)$, solves this optimization problem and
# returns a $u^*$ corresponding to an optimal solution.
#
# !!! note
#     We add $\bar{x}$ as a decision variable, along with the fishing constraint
#     $\bar{x} = x$ for two reasons: it makes it obvious that formulating a
#     problem with $x \times u$ results in a bilinear program instead of a
#     linear program, and it simplifies that internal algorithm that SDDP.jl
#     uses to find an optimal policy.

# These subproblems are very difficult to solve exactly, because they involve
# recursive optimization problems with lots of nested expectations.

# Therefore, instead of solving them exactly, SDDP works by iteratively
# approximating the expectation term of each subproblem, which is also called
# the cost-to-go term. For now, you don't need to understand the details, we
# will explain how shortly.

# The subproblem view of a multistage stochastic program is also important,
# because it provides a convienient way of communicating the different parts of
# the broader problem, and it is how we will communicate the problem to SDDP.jl.
# All we need to do is drop the cost-to-go term and fishing constraint, and
# define a new subproblem `SP` as:
# ```math
# \begin{aligned}
# \texttt{SP}_i(x, \omega) : \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) \\
# & x^\prime = T_i(\bar{x}, u, \omega) \\
# & u \in U_i(\bar{x}, \omega) \\
# \end{aligned}
# ```
# !!! note
#     When we talk about formulating a **subproblem** with SDDP.jl, this is the
#     formulation we mean.

# We've retained the transition function and uncertainty set because they help
# to motivate the different components of the subproblem. However, in general,
# the subproblem can be more general. A better (less restrictive) representation
# might be:
# ```math
# \begin{aligned}
# \texttt{SP}_i(x, \omega) : \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, x^\prime, u, \omega) \\
# & (\bar{x}, x^\prime, u) \in \mathcal{X}_i(\omega)
# \end{aligned}
# ```
# Note that the outgoing state variable can appear in the objective, and we can
# add constraints involving the incoming and outgoing state variables. It
# should be obvious how to map between the two representations.

# ## Preliminaries: Kelley's cutting plane algorithm

# Kelley's cutting plane algorithm is an iterative method for minimizing convex
# functions. Given a convex function $f(x)$, Kelley's constructs an
# under-approximation of the function at the minimum by a set of first-order
# Taylor series approximations (called **cuts**) constructed at a set of $K$
# points $k = 1,\ldots,K$:
# ```math
# \begin{aligned}
# f^K = \min\limits_{\theta \in \mathbb{R}, x \in \mathbb{R}^N} \;\; & \theta\\
# & \theta \ge f(x_k) + \frac{df}{dx}\left(x_k\right)^\top (x - x_k),\quad k=1,\ldots,K\\
# & \theta \ge M,
# \end{aligned}
# ```
# where $M$ is a sufficiently large negative number that is a lower bound for
# $f$ over the domain of $x$.

# As more cuts are added:
# ```math
# \lim_{K \rightarrow \infty} f^K = \min\limits_{x \in \mathbb{R}^N} f(x)
# ```

# ### Bounds

# By convexity, $f^K \le f(x)$ for all $x$. Thus, if $x^*$ is a minimizer of
# $f$, then at any point in time we can construct a lower bound for $f(x^*)$ by
# solving $f^K$.

# Moreover, since any feasible point is an upper bound, we can use the primal
# solution $x^K$ returned by solving $f^K$ to evaluate $f(x_K)$ to generate an
# upper bound.

# Therefore, $f^K \le f(x^*) \le f(x_K)$.

# ### Implementation

# Here is pseudo-code fo the Kelley algorithm:

# 1. Take as input a function $f$ and a iteration limit $K_{max}$. Set $K = 0$,
#    and initialize $f^K$. Set $lb = -\infty$ and $ub = \infty$
# 2. Solve $f^K$ to obtain a candidate solution $x_{K+1}$.
# 3. Update $lb = f^K$ and $ub = \min\{ub, f(x_{K+1}\}$
# 4. Add a cut $\theta \ge f(x_{K+1}) + \frac{df}{dx}\left(x_{K+1}\right)^\top (x - x_{K+1})$ to form $f^{K+1}$.
# 5. Increment $K$
# 6. If $K = K_{max}$, STOP, otherwise, go to step 2.

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
    lower_bound, upper_bound = -Inf, Inf
    while true
        ## Step (2)
        JuMP.optimize!(model)
        x_k = JuMP.value.(x)
        ## Step (3)
        lower_bound = JuMP.objective_value(model)
        upper_bound = min(upper_bound, f(x_k))
        println("K = $K : $(lower_bound) <= f(x*) <= $(upper_bound)")
        ## Step (4):
        c = JuMP.@constraint(model, θ >= f(x_k) + dfdx(x_k)' * (x .- x_k))
        ## Step (5):
        K = K + 1
        ## Step (6):
        if K == iteration_limit
            break
        end
    end
    println("Found solution: x_K = ", JuMP.value.(x))
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

# ## Preliminaries: approximating the cost-to-go term

# In the background theory section, we discussed how you could formulate an
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
# manner. Before we get to that though, let's start writing some code.

# ## Implementation: modeling

# Let's make a start by defining the problem structure. Like SDDP.jl, we need a
# few things:
#
# 1. A description of the structure of the policy graph: how many nodes there
#    are, and the arcs linking the nodes together with their corresponding
#    probabilities.
# 2. A JuMP model for each node in the policy graph
# 3. A way to identify the incoming and outgoing state variables of each node
# 4. A description of the random variable, as well as a function that we can
#    call that will modify the JuMP model to reflect the realization of the
#    random variable.
# 5. A decision variable to act as the approximated cost-to-go term.

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

# Now we have two building blocks, we can declare the structure of each node.

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
# large amount of information to the screen when creating a model.

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

# We're going to copy the example from [Basic II: adding uncertainty](@ref),
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

# Then, we can create a model using the `subproblem_builder` function we defined
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

# ## Implementation: solution algorithm

# Before we get properly coding the solution algorithm, it's also going to be
# useful to have a function that samples a realization of the random variable
# defined by `Ω` and `P`:

function sample_uncertainty(uncertainty::Uncertainty)
    r = rand()
    for (p, ω) in zip(uncertainty.P, uncertainty.Ω)
        r -= p
        if r < 0.0
            return ω
        end
    end
    error("We should never get here because P should sum to 1.0.")
end

# You should be able to work out what is going on. `rand()` samples a uniform
# random variable in `[0, 1)`. For example:

for _ = 1:3
    println("ω = ", sample_uncertainty(model.nodes[1].uncertainty))
end

# It's also going to be useful to define a function that generates a random walk
# through the nodes of the graph:

function sample_next_node(model::PolicyGraph, current::Int)
    r = rand()
    for (to, probability) in model.arcs[current]
        r -= probability
        if r < 0.0
            return to
        end
    end
    ## We looped through the outgoing arcs and still have probability left over!
    ## This means we've hit a leaf node and it's time to stop walking.
    return nothing
end

# For example:

for i = 1:3
    ## We use `repr` to print the next node, because `sample_next_node` can
    ## return `nothing`.
    println("Next node from $(i) = ", repr(sample_next_node(model, i)))
end

# This is a little boring, because our graph is simple. However, more
# complicated graphs will generate more interesting trajectories!

# ## Implementation: the forward pass

# Like Kelley's algorithm, we need a way of generating candidate solutions
# $x_K$. However, unlike the Kelley's example, our functions need two inputs:
# an incoming state variable and a realization of the random variable. We get
# these from a simulation of the policy, which we call the **forward pass**.

# The forward pass walks the policy graph from start to end, transitioning
# randomly along the arcs. At each node, it observes a realization of the random
# variable and solves the approximated subproblem to generate a candidate
# outgoing state variable $x_k^\prime$. The outgoing state variable is passed as
# the incoming state variable to the next node in the trajectory.

function forward_pass(
    model::PolicyGraph,
    io::IO = stdout,
)
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
        stage_cost = JuMP.objective_value(node.subproblem) - JuMP.value(node.cost_to_go)
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
# corresponding outgoing state variables. In the **backward pass**, we walk back
# up this list, and at each node, we compute the cut:
# ```math
# \theta \ge \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi) + \frac{dV_j^k}{dx^\prime}\left(x^\prime_k, \varphi\right)^\top (x^\prime - x^\prime_k)\right],\quad k=1,\ldots,K
# ```

function backward_pass(
    model::PolicyGraph,
    trajectory::Vector{Tuple{Int,Dict{Symbol,Float64}}},
    io::IO = stdout,
)
    println(io, "| Backward pass")
    ## For the backward pass, we walk back up the nodes.
    for i = reverse(1:length(trajectory))
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
        ## Now for each possible node and realization of the uncertainty, solve
        ## the subproblem, and add `P_ij * p_ω * [V + dVdxᵀ(x - x_k)]` to the
        ## cut expression.
        for (j, P_ij) in model.arcs[index]
            next_node = model.nodes[j]
            for (k, v) in outgoing_states
                JuMP.fix(next_node.states[k].in, v; force = true)
            end
            for (pω, ω) in zip(next_node.uncertainty.P, next_node.uncertainty.Ω)
                println(io, "| | | Solving ω = ", ω)
                next_node.uncertainty.parameterize(ω)
                JuMP.optimize!(next_node.subproblem)
                V = JuMP.objective_value(next_node.subproblem)
                println(io, "| | | | V = ", V)
                dVdx = Dict(k => JuMP.reduced_cost(v.in) for (k, v) in next_node.states)
                println(io, "| | | | dVdx′ = ", dVdx)
                cut_expression += P_ij * pω * JuMP.@expression(
                    node.subproblem,
                    V + sum(
                        dVdx[k] * (x.out - outgoing_states[k])
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
    return nothing
end

# ## Implementation: the training loop

# ### Lower bounds

# Recall from Kelley's that we can obtain a lower bound for $f(x^*)$ be
# evaluating $f^K$. The analagous lower bound for a multistage stochastic
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
#    vector.

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
    simulations = [forward_pass(model, devnull) for _ = 1:replications]
    z = [s[2] for s in simulations]
    μ  = Statistics.mean(z)
    tσ = 1.96 * Statistics.std(z) / sqrt(replications)
    return μ, tσ
end

# !!! note
#     The width of the confidence interval is incorrect if there are cycles in
#     the graph, because the distribution of simulation costs `z` is not
#     symmetric. The mean is correct, however.

# ### Training

# The `train` loop of SDDP just applies the forward and backward passes
# iteratively, followed by a final simulation to compute the upper bound
# confidence interval:

function train(
    model::PolicyGraph;
    iteration_limit::Int,
    replications::Int,
    io::IO = stdout,
)
    for i = 1:iteration_limit
        println(io, "Starting iteration $(i)")
        outgoing_states, _ = forward_pass(model, io)
        backward_pass(model, outgoing_states, io)
        println(io, "| Finished iteration")
        println(io, "| | lower_bound = ", lower_bound(model))
    end
    μ, tσ = upper_bound(model; replications = replications)
    println(io, "Upper bound = $(μ) ± $(tσ)")
    return
end

# Using our `model` we defined earlier, we can go:

train(model; iteration_limit = 3, replications = 100)

# Success! We trained a policy for a finite horizon multistage stochastic
# program using stochastic dual dynamic programming.

# ## Implementation: decision rules

# A final step is the ability to evaluate the decision rule associated with a
# node:

function decision_rule(
    node::Node;
    incoming_state::Dict{Symbol,Float64},
    random_variable,
)
    node.uncertainty.parameterize(random_variable)
    for (k, v) in incoming_state
        JuMP.fix(node.states[k].in, v; force = true)
    end
    JuMP.optimize!(node.subproblem)
    return Dict(
        k => JuMP.value.(v)
        for (k, v) in JuMP.object_dictionary(node.subproblem)
    )
end

decision_rule(
    model.nodes[1];
    incoming_state = Dict(:volume => 150.0),
    random_variable = 75,
)

# Note how the random variable can be **out-of-sample**, i.e., it doesn't have
# to be in the vector $\Omega$ we created when defining the model! This is a
# notable difference to other multistage stochastic solution methods like
# progressive hedging or using the deterministic equivalent.

# ## Example: infinite horizon

# As promised earlier, our implementation is actually pretty general. It can
# solve any multistage stochastic (linear) program defined by a policy graph,
# including infinite horizon problems!

# Here's an example, where we have extended our earlier problem with an arc from
# node 3 to node 2 with probability 0.5. You can interpret the 0.5 as a discount
# factor.

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

decision_rule(
    model.nodes[3];
    incoming_state = Dict(:volume => 100.0),
    random_variable = 10.0,
)
