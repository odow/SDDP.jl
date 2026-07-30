# # Hydro-thermal scheduling

# ## Problem Description

# In a hydro-thermal problem, the agent controls a hydro-electric generator and reservoir.
# Each time period, they need to choose a generation quantity from thermal `g_t`, and hydro
# `g_h`, in order to meed demand `w_d`, which is a stagewise-independent random variable.
# The state variable, `x`, is the quantity of water in the reservoir at the start of each
# time period, and it has a minimum level of 5 units and a maximum level of 15 units. We
# assume that there are 10 units of water in the reservoir at the start of time, so that
# `x_0 = 10`. The state-variable is connected through time by the water balance constraint:
# `x.out = x.in - g_t - s + w_i,` where `x.out` is the quantity of water at the end of the
# time period, `x.in` is the quantity of water at the start of the time period, `s` is the
# quantity of water spilled from the reservoir, and `w_i` is a stagewise-independent random
# variable that represents the inflow into the reservoir during the time period.

# We assume that there are three stages, `t=1, 2, 3`, representing summer-fall, winter, and
# spring, and that we are solving this problem in an infinite-horizon setting with a
# discount factor of `0.95`.

# In each stage, the agent incurs the cost of spillage, plus the cost of thermal generation.
# We assume that the cost of thermal generation is dependent on the stage `t = 1, 2, 3`, and
# that in each stage, `w` is drawn from the set `(w_i, w_d) = {(0, 7.5), (3, 5), (10, 2.5)}`
# with equal probability.

# ## Importing packages

# For this example, in addition to `SDDP`, we need `GLPK` as a solver and `Statisitics` to
# compute the mean of our simulations.

using GLPK
using SDDP
using Statistics

# ## Constructing the policy graph

# There are three stages in our problem, so we construct a linear policy graph with three
# stages using [`SDDP.LinearGraph`](@ref):

graph = SDDP.LinearGraph(3)

# Then, because we want to solve an infinite-horizon problem, we add an additional edge
# between node `3` and node `1` with probability `0.95`:

SDDP.add_edge(graph, 3 => 1, 0.95)

# ## Constructing the model

# Much of the macro code (i.e., lines starting with `@`) in the first part of the following
# should be familiar to users of JuMP.

# Inside the `do-end` block (if this isn't familiar, see
# [What's this weird `do` syntax?](@ref)), `sp` is a standard JuMP model, and `t` is an
# index for the state variable that will be called with `t = 1, 2, 3`.

# The state variable `x`, constructed by passing the `SDDP.State` tag to `@variable` is
# actually a Julia struct with two fields: `x.in` and `x.out` corresponding to the incoming
# and outgoing state variables respectively. Both `x.in` and `x.out` are standard JuMP
# variables. The `initial_value` keyword provides the value of the state variable in the
# root node (i.e., `x_0`).

# Compared to a JuMP model, one key difference is that we use [`@stageobjective`](@ref)
# instead of `@objective`. The [`SDDP.parameterize`](@ref) function takes a list of supports
# for `w` and parameterizes the JuMP model `sp` by setting the right-hand sides of the
# appropriate constraints (note how the constraints initially have a right-hand side of
# `0`). By default, it is assumed that the realizations have uniform probability, but a
# probability mass vector can also be provided.

model = SDDP.PolicyGraph(
    graph, sense = :Min, lower_bound = 0.0, optimizer = GLPK.Optimizer
) do sp, t
    @variable(sp, 5 <= x <= 15, SDDP.State, initial_value = 10)
    @variable(sp, g_t >= 0)
    @variable(sp, g_h >= 0)
    @variable(sp, s >= 0)
    @constraint(sp, balance, x.out - x.in + g_h + s == 0)
    @constraint(sp, demand, g_h + g_t == 0)
    @stageobjective(sp, s + t * g_t)
    SDDP.parameterize(sp, [[0, 7.5], [3, 5], [10, 2.5]]) do w
        set_normalized_rhs(balance, w[1])
        set_normalized_rhs(demand, w[2])
    end
end

# ## Training the policy

# Once a model has been constructed, the next step is to train the policy. This can be
# achieved using [`SDDP.train`](@ref). There are many options that can be passed, but
# `iteration_limit` terminates the training after the prescribed number of SDDP iterations.

SDDP.train(model, iteration_limit = 100)

# ## Simulating the policy

# After training, we can simulate the policy using [`SDDP.simulate`](@ref).

sims = SDDP.simulate(model, 100, [:g_t])
mu = round(mean([s[1][:g_t] for s in sims]), digits = 2)
println("On average, $(mu) units of thermal are used in the first stage.")

# ## Extracting the water values

# Finally, we can use [`SDDP.ValueFunction`](@ref) and [`SDDP.evaluate`](@ref) to obtain and
# evaluate the value function at different points in the state-space. Note that since we
# are minimizing, the price has a negative sign: each additional unit of water leads to a
# decrease in the the expected long-run cost.

V = SDDP.ValueFunction(model[1])
cost, price = SDDP.evaluate(V, x = 10)
