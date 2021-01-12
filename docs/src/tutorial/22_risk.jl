# # Theory II: risk aversion

# In [Theory I: an intro to SDDP](@ref), we implemented a basic version of the
# SDDP algorithm.

# This tutorial extends that implementation to add **risk-aversion**.

# ## Risk aversion: what and why?

# Often, the agents making decisions in complex systems are **risk-averse**,
# that is, they care more about avoiding very bad outcomes, than they do about
# having a good average outcome.

# As an example, consumers in a hydro-thermal problem may be willing to pay a
# slightly higher electricity price on average, if it means that there is a
# lower probability of blackouts.

# Risk aversion in multistage stochastic programming has been well studied in
# the academic literature, and is widely used in production implementations
# around the world.

# ## Risk measures

# One way to add risk aversion to models is to use a **risk measure**. A risk
# measure is a function that maps a random variable to a real number.
#
# You are probably already familiar with lots of different risk measures. For
# example, the mean, median, mode, and maximum are all risk measures.
#
# We call the act of applying a risk measure to a random variable "computing the
# risk" of a random variable.
#
# To keep things simple, and because we need it for SDDP, we restrict our
# attention to random variables $Z$ with a finite sample space $\Omega$
# and positive probabilities $p_\omega$ for all $\omega \in \Omega$. We denote
# the realizations of $Z$ by $Z(\omega) = z_\omega$.

# A risk measure, $\mathbb{F}[Z]$, is a **convex risk measure** if it satisfies
# the following axioms:
#
# **Axiom 1: monotonicity**
#
# Given two random variables $Z_1$ and $Z_2$, with $Z_1 \le Z_2$ almost surely,
# then $\mathbb{F}[Z_1] \le F[Z_2]$.
#
# **Axiom 2: translation equivariance**
#
# Given two random variables $Z_1$ and $Z_2$, then for all $a \in \mathbb{R}$,
# $\mathbb{F}[Z + a] = \mathbb{F}[Z] + a$.
#
# **Axiom 3: convexity**
#
# Given two random variables $Z_1$ and $Z_2$, then for all $a \in [0, 1]$,
# $\mathbb{F}[a Z_1 + (1 - a) Z_2] \le a \mathbb{F}[Z_1] + (1-a)\mathbb{F}[Z_2]$.

# Now we know what a risk measure is, let's see how we can use them to form
# risk-averse decision rules.

# ## Risk-averse decision rules: Part I

# We started this tutorial by explaining that we are interested in risk aversion
# because some agents are risk-averse. What that really means, is that they
# want a policy that is also risk-averse. The question then becomes, how do we
# create risk-averse decision rules and policies?

# Recall from [Theory I: an intro to SDDP](@ref) that we can form an optimal
# decision rule using the recursive formulation:
# ```math
# \begin{aligned}
# V_i(x, \omega) = \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) + \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}[V_j(x^\prime, \varphi)]\\
# & x^\prime = T_i(\bar{x}, u, \omega) \\
# & u \in U_i(\bar{x}, \omega) \\
# & \bar{x} = x,
# \end{aligned}
# ```
# where our decision rule, $\pi_i(x, \omega)$, solves this optimization problem
# and returns a $u^*$ corresponding to an optimal solution.

# If we can replace the expectation operator $\mathbb{E}$ with another (more
# risk-averse) risk measure $\mathbb{F}$, then our decision rule will attempt to
# choose a control decision now that minimizes the risk of the future costs, as
# opposed to the expectation of the future costs.

# ## Primal risk measures

# Now we know what a risk measure is, and how we will use it, let's implement
# some code to see how we can compute the risk of some random variables.

# For reasons that will become apparent shortly, we will implement the primal
# version of each risk measure. (We will discuss the dual version in the next
# section.)

# First, we need some data:

Z = [1.0, 2.0, 3.0, 4.0]

# with probabilities:

p = [0.1, 0.2, 0.4, 0.3]

# We're going to implement a number of different risk measures, so to leverage
# Julia's multiple dispatch, we create an abstract type:

abstract type AbstractRiskMeasure end

# and function to overload:

"""
    primal_risk(F::AbstractRiskMeasure, Z::Vector{Float64}, p::Vector{Float64})

Use `F` to compute the risk of the random variable defined by a vector of costs
`Z` and non-zero probabilities `p`.
"""
function primal_risk end

# ### Expectation

# The expectation, $\mathbb{E}$, also called the mean or the average, is the
# most widely used convex risk measure. The expectation of a random variable is
# just the sum of $Z$ weighted by the probability:

struct Expectation <: AbstractRiskMeasure end

function primal_risk(::Expectation, Z::Vector{Float64}, p::Vector{Float64})
    return sum(p[i] * Z[i] for i = 1:length(p))
end

# Let's try it out:

primal_risk(Expectation(), Z, p)

# ### WorstCase

# The worst-case risk measure, also called the maximum, is another widely used
# convex risk measure. This risk measure doesn't care about the probability
# vector `p`, only the cost vector `Z`:

struct WorstCase <: AbstractRiskMeasure end

function primal_risk(::WorstCase, Z::Vector{Float64}, ::Vector{Float64})
    return maximum(Z)
end

# Let's try it out:

primal_risk(WorstCase(), Z, p)

# ### Entropic

# A more interesting, and less widely used risk measure is the entropic risk
# measure. The entropic risk measure is parameterized by a value $\gamma > 0$,
# and computes the risk of a random variable as:
# ```math
# \mathbb{F}_\gamma[Z] = \frac{1}{\gamma}\log\left(\mathbb{E}[e^{\gamma Z}]\right).
# ```

struct Entropic <: AbstractRiskMeasure
    γ::Float64
    function Entropic(γ)
        if !(γ > 0)
            throw(DomainError(γ, "Entropic risk measure must have γ > 0."))
        end
        return new(γ)
    end
end

function primal_risk(F::Entropic, Z::Vector{Float64}, p::Vector{Float64})
    FZ = 1 / F.γ * log(sum(p[i] * exp(F.γ * big(Z[i])) for i = 1:length(p)))
    return Float64(FZ)
end

# !!! warning
#     We use `big(Z[i])` in the `exp` call because `exp(x)` overflows when
#     $x > 709$. We explicitly cast back to `Float64` for the `return` value.

# Let's try it out for different values of $\gamma$:

for γ in [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1_000.0]
    FZ = primal_risk(Entropic(γ), Z, p)
    println("γ = $(γ), F[Z] = $(FZ)")
end

# !!! info
#     The entropic has two extremes. As $\gamma \rightarrow 0$, the entropic
#     acts like the expectation risk measure, and as $\gamma \rightarrow \infty$,
#     the entropic acts like the worst-case risk measure.

# Computing risk measures this way work for fixed data, but it's not so apparent
# how we can implement this into a linear program, because our cut
# calculation requires:
# ```math
# \theta \ge \mathbb{F}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi)\right] + \frac{d}{dx^\prime}\mathbb{F}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi)\right]^\top (x^\prime - x^\prime_k).
# ```

# !!! warning
#     Note how we need to expliclty compute the risk-averse gradient! When
#     constucting the cuts with the expectation operator, we implicitly used the
#     law of total expectation to combine the two expectations; we can't do that
#     for a general risk measure.

# There is a nice solution to this problem, and that is to use the dual
# representation of a risk measure, instead of the primal.

# ## Dual risk measures

# Convex risk measures have a dual representation as follows:
# ```math
# \mathbb{F}[Z] = \sup\limits_{q \in\mathcal{M}(p)} \mathbb{E}_q[Z] - \alpha(p, q),
# ```
# where $\mathcal{M}(p)$ is a convex set, and $\alpha$ is a function that maps
# the probability vectors $p$ and $q$ to a real number.

# The dual of a convex risk measure can be interpreted as taking the expectation
# of the random variable $Z$ with respect to the worst probability vector $q$
# that lies within the set $\mathcal{M}$, less some concave penalty term
# $\alpha(p, q)$.

# If we define a function `dual_risk_inner` that computes `q` and `α`:

"""
    dual_risk_inner(
        F::AbstractRiskMeasure, Z::Vector{Float64}, p::Vector{Float64}
    )::Tuple{Vector{Float64},Float64}

Return a tuple formed by the worst-case probability vector `q` and the
corresponding evaluation `α(p, q)`.
"""
function dual_risk_inner end

# then we can write a generic `dual_risk` function as:

function dual_risk(
    F::AbstractRiskMeasure,
    Z::Vector{Float64},
    p::Vector{Float64},
)
    q, α = dual_risk_inner(F, Z, p)
    return sum(q[i] * Z[i] for i = 1:length(q)) - α
end

# ### Expectation

# For expectation, $\mathcal{M}(p) = \{p\}$, and $\alpha(\cdot, \cdot) = 0$.
# Therefore:

function dual_risk_inner(::Expectation, ::Vector{Float64}, p::Vector{Float64})
    return p, 0.0
end

# We can check we get the same result as the primal version:

F = Expectation()
dual_risk(F, Z, p) == primal_risk(F, Z, p)

# ### Worst-case

# For worst-case, $\mathcal{M}(p) = \{q > 0\;|\;\sum\limits_{\omega\in\Omega}q_\omega = 1\}$,
# and $\alpha(\cdot, \cdot) = 0$. Therefore, the dual representation just puts
# all of the probability weight on the maximum outcome:

function dual_risk_inner(::WorstCase, Z::Vector{Float64}, ::Vector{Float64})
    q = zeros(length(Z))
    _, index = findmax(Z)
    q[index] = 1.0
    return q, 0.0
end

# We can check we get the same result as the primal version:

F = WorstCase()
dual_risk(F, Z, p) == primal_risk(F, Z, p)

# ### Entropic

# Computing the $q^*$ associated with the entropic risk measure is a little
# more complicated. However, [Dowson, Morton, and Pagnoncelli (2020)]() show
# that:
# ```math
# q_\omega^* = \frac{p_\omega e^{\gamma z_\omega}}{\sum\limits_{\varphi \in \Omega} p_\varphi e^{\gamma z_\varphi}}.
# ```
# In addition:
# ```math
# \alpha(p, q) = \sum\limits_{\omega\in\Omega} q_\omega \log\left(\frac{q_\omega}{p_\omega}\right).
# ```

function dual_risk_inner(F::Entropic, Z::Vector{Float64}, p::Vector{Float64})
    q, α = zeros(length(p)), big(0.0)
    peγz = p .* exp.(F.γ .* big.(Z))
    sum_peγz = sum(peγz)
    for i = 1:length(q)
        big_q = peγz[i] / sum_peγz
        α += big_q * log(big_q / p[i])
        q[i] = Float64(big_q)
    end
    return q, Float64(α / F.γ)
end

# !!! warning
#     Again, note that we use `big` to avoid introducing overflow errors, before
#     explicitly casting back to `Float64` for the values we return.

# We can check we get the same result as the primal version:

for γ in [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1_000.0]
    primal = primal_risk(Entropic(γ), Z, p)
    dual = dual_risk(Entropic(γ), Z, p)
    success = primal ≈ dual ? "✓" : "×"
    println("$(success) γ = $(γ), primal = $(primal), dual = $(dual)")
end

# ### Risk-averse gradients

# We ended the section on primal risk measures by explaining how we couldn't
# use the primal risk measure in the cut calculation because we needed some way
# of computing the risk-averse gradient:
# ```math
# \theta \ge \mathbb{F}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi)\right] + \frac{d}{dx^\prime}\mathbb{F}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi)\right]^\top (x^\prime - x^\prime_k),
# ```

# The reason we use the dual representation is because of the following theorem,
# which explains how to compute the risk-averse gradient:

# !!! info "Theorem"
#     Let $\omega \in \Omega$ index a random vector with finite support and with
#     nominal probability mass function, $p \in \mathcal{P}$, which satisfies
#     $p > 0$.
#
#     Consider a convex risk measure, $\mathbb{F}$, with a convex risk set,
#     $\mathcal{M}(p)$, so that $\mathbb{F}$ can be expressed as the dual form.
#
#     Let $V(x,\omega)$ be convex with respect to $x$ for all fixed
#     $\omega\in\Omega$, and let $\lambda(\tilde{x}, \omega)$ be a subgradient
#     of $V(x,\omega)$ with respect to $x$ at $x = \tilde{x}$ for each
#     $\omega \in \Omega$.
#
#     Then, $\sum_{\omega\in\Omega}q^*_{\omega} \lambda(\tilde{x},\omega)$ is a
#     subgradient of $\mathbb{F}[V(x,\omega)]$ at $\tilde{x}$, where
#     $q^* \in \argmax_{q \in \mathcal{M}(p)}\left\{{\mathbb{E}}_q[V(\tilde{x},\omega)] - \alpha(p, q)\right\}$.

# This theorem can be a little hard to unpack, so let's see an example:

import ForwardDiff

function risk_averse_subgradient(
    V::Function,
    ## Use automatic differentiation to compute the gradient of V w.r.t. x,
    ## given a fixed ω.
    λ::Function = (x, ω) -> ForwardDiff.gradient(x -> V(x, ω), x);
    F::AbstractRiskMeasure,
    Ω::Vector,
    p::Vector{Float64},
    x̃::Vector{Float64},
)
    V_ω = [V(x̃, ω) for ω in Ω]
    q, α = dual_risk_inner(F, V_ω, p)
    dVdx = sum(q[i] * λ(x̃, ω) for (i, ω) in enumerate(Ω))
    return q' * V_ω - α, dVdx
end

# As our example, we use

V(x, ω) = ω * x[1]^2

# with

Ω = [1.0, 2.0, 3.0]

#  and

p = [0.3, 0.4, 0.3]

# at the point

x̃ = [3.0]

# If $\mathbb{F}$ is the expectation operator, then
# ```math
# \mathbb{E}[V(x, \omega)] =  2 x^2
# ```
# and so the function evaluation $x=3$ is $18$ and the subgradient is $12$.
# Let's check we get it right:

risk_averse_subgradient(V; F = Expectation(), Ω = Ω, p = p, x̃ = x̃)

# If $\mathbb{F}$ is the worst-case, then
# ```math
# \mathbb{F}[V(x, \omega)] = 3 x^2
# ```
# and so the function evaluation at $x=3$ is $27$, and the subgradient
# is $18$. Let's check we get it right:

risk_averse_subgradient(V; F = WorstCase(), Ω = Ω, p = p, x̃ = x̃)

# If $\mathbb{F}$ is the entropic, the math is a little more difficult,
# so you'll have to trust that the math is correct:

for γ in [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1_000.0]
    ∇V = risk_averse_subgradient(V; F = Entropic(γ), Ω = Ω, p = p, x̃ = x̃)
    println("γ = $(γ), dF[V(x, ω)]/dx  = $(∇V)")
end

# For a sanity check, as $\gamma \rightarrow 0$, we get the subgradient of the
# expectation operator, and as $\gamma \rightarrow \infty$, we get the
# subgradient of the worse-case risk measure. You can verify the values
# analytically using the primal definition of the entropic risk measure.

# # Risk-averse decision rules: Part II

# Using the dual representation of a convex risk measure, we can re-write the
# cut:
# ```math
# \theta \ge \mathbb{F}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi)\right] + \frac{d}{dx^\prime}\mathbb{F}_{j \in i^+, \varphi \in \Omega_j}\left[V_j^k(x^\prime_k, \varphi)\right]^\top (x^\prime - x^\prime_k),\quad k=1,\ldots,K
# ```
# as:
# ```math
# \theta \ge \mathbb{E}_{q_k}\left[V_j^k(x^\prime_k, \varphi) + \frac{d}{dx^\prime}V_j^k(x^\prime_k, \varphi)^\top (x^\prime - x^\prime_k)\right] - \alpha(p, q_k),\quad k=1,\ldots,K,
# ```
# where $q_k = \mathrm{arg}\sup\limits_{q \in\mathcal{M}(p)} \mathbb{E}_q[V_j^k(x_k^\prime, \varphi)] - \alpha(p, q)$

# Therefore, we can formulate a risk-averse decision rule as:
# ```math
# \begin{aligned}
# V_i^K(x, \omega) = \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) + \theta\\
# & x^\prime = T_i(\bar{x}, u, \omega) \\
# & u \in U_i(\bar{x}, \omega) \\
# & \bar{x} = x \\
# & \theta \ge \mathbb{E}_{q_k}\left[V_j^k(x^\prime_k, \varphi) + \frac{d}{dx^\prime}V_j^k(x^\prime_k, \varphi)^\top (x^\prime - x^\prime_k)\right] - \alpha(p, q_k),\quad k=1,\ldots,K \\
# & \theta \ge M.
# \end{aligned}
# ```
# where $q_k = \mathrm{arg}\sup\limits_{q \in\mathcal{M}(p)} \mathbb{E}_q[V_j^k(x_k^\prime, \varphi)] - \alpha(p, q)$

# ## Implementation

# Now we're ready to implement our risk-averse version of SDDP.

# As a prerequisite, we need most of the code from [Theory I: an introduction to SDDP](@ref).

# ```@raw html
# <details>
# <summary>Click to view code from the tutorial "Theory I: an introduction to SDDP".</summary>
# ```

import GLPK
import JuMP
import Statistics

struct State
    in::JuMP.VariableRef
    out::JuMP.VariableRef
end

struct Uncertainty
    parameterize::Function
    Ω::Vector{Any}
    P::Vector{Float64}
end

struct Node
    subproblem::JuMP.Model
    states::Dict{Symbol,State}
    uncertainty::Uncertainty
    cost_to_go::JuMP.VariableRef
end

struct PolicyGraph
    nodes::Vector{Node}
    arcs::Vector{Dict{Int,Float64}}
end

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

function PolicyGraph(
    subproblem_builder::Function;
    graph::Vector{Dict{Int,Float64}},
    lower_bound::Float64,
    optimizer,
)
    nodes = Node[]
    for t = 1:length(graph)
        model = JuMP.Model(optimizer)
        states, uncertainty = subproblem_builder(model, t)
        JuMP.@variable(model, cost_to_go >= lower_bound)
        obj = JuMP.objective_function(model)
        JuMP.@objective(model, Min, obj + cost_to_go)
        if length(graph[t]) == 0
            JuMP.fix(cost_to_go, 0.0; force = true)
        end
        push!(nodes, Node(model, states, uncertainty, cost_to_go))
    end
    return PolicyGraph(nodes, graph)
end

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

function sample_next_node(model::PolicyGraph, current::Int)
    if length(model.arcs[current]) == 0
        return nothing
    else
        r = rand()
        for (to, probability) in model.arcs[current]
            r -= probability
            if r < 0.0
                return to
            end
        end
        return nothing
    end
end

function forward_pass(model::PolicyGraph, io::IO = stdout)
    incoming_state = Dict(
        k => JuMP.fix_value(v.in) for (k, v) in model.nodes[1].states
    )
    simulation_cost = 0.0
    trajectory = Tuple{Int,Dict{Symbol,Float64}}[]
    t = 1
    while t !== nothing
        node = model.nodes[t]
        ω = sample_uncertainty(node.uncertainty)
        node.uncertainty.parameterize(ω)
        for (k, v) in incoming_state
            JuMP.fix(node.states[k].in, v; force = true)
        end
        JuMP.optimize!(node.subproblem)
        if JuMP.termination_status(node.subproblem) != JuMP.MOI.OPTIMAL
            error("Something went terribly wrong!")
        end
        outgoing_state = Dict(k => JuMP.value(v.out) for (k, v) in node.states)
        stage_cost = JuMP.objective_value(node.subproblem) - JuMP.value(node.cost_to_go)
        simulation_cost += stage_cost
        incoming_state = outgoing_state
        push!(trajectory, (t, outgoing_state))
        t = sample_next_node(model, t)
    end
    return trajectory, simulation_cost
end

function upper_bound(model::PolicyGraph; replications::Int)
    simulations = [forward_pass(model, devnull) for i = 1:replications]
    z = [s[2] for s in simulations]
    μ  = Statistics.mean(z)
    tσ = 1.96 * Statistics.std(z) / sqrt(replications)
    return μ, tσ
end

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
        k => JuMP.value.(v)
        for (k, v) in JuMP.object_dictionary(the_node.subproblem)
    )
end

# ```@raw html
# </details>
# ```

# We only need to change two functions; `backward_pass` and `train`.

# And here's the backward pass

function backward_pass(
    model::PolicyGraph,
    trajectory::Vector{Tuple{Int,Dict{Symbol,Float64}}},
    io::IO = stdout;
    risk_measure::AbstractRiskMeasure,
)
    println(io, "| Backward pass")
    for i = reverse(1:length(trajectory))
        index, outgoing_states = trajectory[i]
        node = model.nodes[index]
        println(io, "| | Visiting node $(index)")
        if length(model.arcs[index]) == 0
            continue
        end
        ## =====================================================================
        ## New! Create vectors to store the cut expressions, V(x,ω) and p:
        cut_expressions, V_ω, p = JuMP.AffExpr[], Float64[], Float64[]
        ## =====================================================================
        for (j, P_ij) in model.arcs[index]
            next_node = model.nodes[j]
            for (k, v) in outgoing_states
                JuMP.fix(next_node.states[k].in, v; force = true)
            end
            for (pφ, φ) in zip(next_node.uncertainty.P, next_node.uncertainty.Ω)
                next_node.uncertainty.parameterize(φ)
                JuMP.optimize!(next_node.subproblem)
                V = JuMP.objective_value(next_node.subproblem)
                dVdx = Dict(
                    k => JuMP.reduced_cost(v.in) for (k, v) in next_node.states
                )
                ## =============================================================
                ## New! Instead of building a single cut expression, append the
                ## expression to cut_expressions:
                push!(
                    cut_expressions,
                    JuMP.@expression(
                        node.subproblem,
                        V + sum(
                            dVdx[k] * (x.out - outgoing_states[k])
                            for (k, x) in node.states
                        ),
                    )
                )
                ## the objective value to Z:
                push!(V_ω, V)
                ## and the probability to p:
                push!(p, P_ij * pφ)
                ## =============================================================
            end
        end
        ## =====================================================================
        ## New! Using the solutions in V_ω, compute q and α:
        q, α = dual_risk_inner(risk_measure, V_ω, p)
        println(io, "| | | Z = ", Z)
        println(io, "| | | p = ", p)
        println(io, "| | | q = ", q)
        println(io, "| | | α = ", α)
        ## Then add the cut:
        c = JuMP.@constraint(
            node.subproblem,
            node.cost_to_go >= sum(
                q[i] * cut_expressions[i] for i = 1:length(q)
            ) - α
        )
        ## =====================================================================
        println(io, "| | | Adding cut : ", c)
    end
    return nothing
end

# We also need to update the `train` loop of SDDP to pass a risk measure to the
# backward pass:

function train(
    model::PolicyGraph;
    iteration_limit::Int,
    replications::Int,
    risk_measure::AbstractRiskMeasure,
    io::IO = stdout,
)
    for i = 1:iteration_limit
        println(io, "Starting iteration $(i)")
        outgoing_states, _ = forward_pass(model, io)
        backward_pass(model, outgoing_states, io; risk_measure = risk_measure)
        println(io, "| Finished iteration")
        println(io, "| | lower_bound = ", lower_bound(model))
    end
    μ, tσ = upper_bound(model; replications = replications)
    println(io, "Upper bound = $(μ) ± $(tσ)")
    return
end

# ## Example: risk-averse hydro-thermal scheduling

# Create the model:

model = PolicyGraph(
    graph = [
        Dict(2 => 1.0),
        Dict(3 => 1.0),
        Dict{Int,Float64}(),
    ],
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
) do subproblem, t
    JuMP.@variable(subproblem, volume_in == 200)
    JuMP.@variable(subproblem, 0 <= volume_out <= 200)
    states = Dict(:volume => State(volume_in, volume_out))
    JuMP.@variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        inflow
    end)
    JuMP.@constraints(subproblem, begin
        volume_out == volume_in + inflow - hydro_generation - hydro_spill
        demand_constraint, thermal_generation + hydro_generation == 150.0
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    JuMP.@objective(subproblem, Min, fuel_cost[t] * thermal_generation)
    uncertainty = Uncertainty([0.0, 50.0, 100.0], [1 / 3, 1 / 3, 1 / 3]) do ω
        JuMP.fix(inflow, ω)
    end
    return states, uncertainty
end

# Then we train a risk-averse policy, passing a risk measure to `train`:

train(
    model;
    iteration_limit = 3,
    replications = 100,
    risk_measure = Entropic(1.0),
)

# Finally, evaluate the decision rule:

evaluate_policy(
    model;
    node = 1,
    incoming_state = Dict(:volume => 150.0),
    random_variable = 75,
)
