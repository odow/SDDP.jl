# # Auto-regressive stochastic processes

# SDDP.jl assumes that the random variable in each node is independent of the
# random variables in all other nodes. However, a common request is to model
# the random variables by some auto-regressive process.

# There are two ways to do this:
#  1. model the random variable as a Markov chain
#  2. use the "state-space expansion" trick

# !!! info
#     This tutorial is in the context of a hydro-thermal scheduling example, but
#     it should be apparent how the ideas transfer to other applications.

using SDDP
import GLPK

# ## [The state-space expansion trick](@id state-space-expansion)

# In [An introduction to SDDP.jl](@ref), we assumed that the inflows were
# stagewise-independent. However, in many cases this is not correct, and inflow
# models are more accurately described by an autoregressive process such as:
# ```math
# inflow_{t} = inflow_{t-1} + \varepsilon
# ```
# Here ``\varepsilon`` is a random variable, and the inflow in stage ``t`` is
# the inflow in stage ``t-1`` plus ``\varepsilon`` (which might be negative).

# For simplicitly, we omit any coefficients and other terms, but this could
# easily be extended to a model like
# ```math
# inflow_{t} = \alpha inflow_{t-1} + \beta + \varepsilon
# ```

# In practice, you can estimate a distribution for ``\varepsilon`` by fitting
# the chosen statistical model to historical data, and then using the empirical
# residuals.

# To implement the auto-regressive model in SDDP.jl, we introduce `inflow` as a
# state variable.

# !!! tip
#     Our rule of thumb for "when is something a state variable?" is: if you
#     need the value of a variable from a previous stage to compute something in
#     stage ``t``, then that variable is a state variable.

model = SDDP.LinearPolicyGraph(
    stages = 3,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
) do sp, t
    @variable(sp, 0 <= x <= 200, SDDP.State, initial_value = 200)
    @variable(sp, g_t >= 0)
    @variable(sp, g_h >= 0)
    @variable(sp, s >= 0)
    @constraint(sp, g_h + g_t == 150)
    c = [50, 100, 150]
    @stageobjective(sp, c[t] * g_t)
    ## =========================================================================
    ## New stuff below Here
    ## Add inflow as a state
    @variable(sp, inflow, SDDP.State, initial_value = 50.0)
    ## Add the random variable as a control variable
    @variable(sp, ε)
    ## The equation describing our statistical model
    @constraint(sp, inflow.out == inflow.in + ε)
    ## The new water balance constraint using the state variable
    @constraint(sp, x.out == x.in - g_h - s + inflow.out)
    ## Assume we have some empirical residuals:
    Ω = [-10.0, 0.1, 9.6]
    SDDP.parameterize(sp, Ω) do ω
        return JuMP.fix(ε, ω)
    end
end

# ### When can this trick be used?

# The state-space expansion trick should be used when:
#
#  * The random variable appears additively in the objective or in the
#    constraints. Something like `inflow * decision_variable` will _not_ work.
#  * The statistical model is linear, or can be written using the JuMP
#    `@constraint` macro.
#  * The dimension of the random variable is small (see
#    [Vector auto-regressive models](@ref) for the multi-variate case).

# ## The Markov chain approach

# In the Markov chain approach, we model the stochastic process for inflow by a
# discrete Markov chain. Markov chains are nodes with transition probabilities
# between the nodes. SDDP.jl has good support for solving problems in which the
# uncertainty is formulated as a Markov chain.

# The first step of the Markov chain approach is to write a function which
# simulates the stochastic process. Here is a simulator for our inflow model:

function simulator()
    inflow = zeros(3)
    current = 50.0
    Ω = [-10.0, 0.1, 9.6]
    for t in 1:3
        current += rand(Ω)
        inflow[t] = current
    end
    return inflow
end

# When called with no arguments, it produces a vector of inflows:

simulator()

# !!! warning
#     The `simulator` must return a `Vector{Float64}`, so it is limited to a
#     uni-variate random variable. It is possible to do something similar for
#     multi-variate random variable, but you'll have to manually construct the
#     Markov transition matrix, and solution times scale poorly, even in the
#     two-dimensional case.

# The next step is to call [`SDDP.MarkovianGraph`](@ref) with our simulator.
# This function will attempt to fit a Markov chain to the stochastic process
# produced by your `simulator`. There are two key arguments:
#  * `budget` is the total number of nodes we want in the Markov chain
#  * `scenarios` is a limit on the number of times we can call `simulator`

graph = SDDP.MarkovianGraph(simulator; budget = 8, scenarios = 30)

# Here we can see we have created a MarkovianGraph with nodes like `(2, 59.7)`.
# The first element of each node is the stage, and the second element is the
# inflow.

# Create a [`SDDP.PolicyGraph`](@ref) using `graph` as follows:

model = SDDP.PolicyGraph(
    graph,  # <--- New stuff
    sense = :Min,
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
) do sp, node
    t, inflow = node  # <--- New stuff
    @variable(sp, 0 <= x <= 200, SDDP.State, initial_value = 200)
    @variable(sp, g_t >= 0)
    @variable(sp, g_h >= 0)
    @variable(sp, s >= 0)
    @constraint(sp, g_h + g_t == 150)
    c = [50, 100, 150]
    @stageobjective(sp, c[t] * g_t)
    ## The new water balance constraint using the node:
    @constraint(sp, x.out == x.in - g_h - s + inflow)
end

# ## When can this trick be used?

# The Markov chain approach should be used when:
#
#  * The random variable is uni-variate
#  * The random variable appears in the objective function or as a variable
#    coefficient in the constraint matrix
#  * It's non-trivial to write the stochastic process as a series of constraints
#    (for example, it uses nonlinear terms)
#  * The number of nodes is modest (for example, a budget of hundreds, up to
#    perhaps 1000)

# ## Vector auto-regressive models

# The [state-space expansion](@ref state-space-expansion) section assumed that
# the random variable was uni-variate. However, the approach naturally extends
# to vector auto-regressive models. For example, if `inflow` is a 2-dimensional
# vector, then we can model a vector auto-regressive model to it as follows:
# ```math
# inflow_{t} = A inflow_{t-1} + b + \varepsilon
# ```
# Here `A` is a 2-by-2 matrix, and `b` and ``\varepsilon`` are 2-by-1 vectors.

 model = SDDP.LinearPolicyGraph(
    stages = 3,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
) do sp, t
    @variable(sp, 0 <= x <= 200, SDDP.State, initial_value = 200)
    @variable(sp, g_t >= 0)
    @variable(sp, g_h >= 0)
    @variable(sp, s >= 0)
    @constraint(sp, g_h + g_t == 150)
    c = [50, 100, 150]
    @stageobjective(sp, c[t] * g_t)
    ## =========================================================================
    ## New stuff below Here
    ## Add inflow as a state
    @variable(sp, inflow[1:2], SDDP.State, initial_value = 50.0)
    ## Add the random variable as a control variable
    @variable(sp, ε[1:2])
    ## The equation describing our statistical model
    A = [0.8 0.2; 0.2 0.8]
    inflow_in = [inflow[i].in for i in 1:2]
    inflow_out = [inflow[i].in for i in 1:2]
    @constraint(sp, inflow_out .== A * inflow_in .+ ε)
    ## The new water balance constraint using the state variable
    @constraint(sp, x.out == x.in - g_h - s + inflow[1].out + inflow[2].out)
    ## Assume we have some empirical residuals:
    Ω₁ = [-10.0, 0.1, 9.6]
    Ω₂ = [-10.0, 0.1, 9.6]
    Ω = [(ω₁, ω₂) for ω₁ in Ω₁ for ω₂ in Ω₂]
    SDDP.parameterize(sp, Ω) do ω
        JuMP.fix(ε[1], ω[1])
        JuMP.fix(ε[2], ω[2])
        return
    end
end
