```@meta
CurrentModule = Kokako
```

# Intermediate IV: objective states

There are many applications in which we want to model a price process that
follows some auto-regressive process. Common examples include stock prices on
financial exchanges and spot-prices in energy markets.

However, it is well known that these cannot be incorporated in to SDDP because
they result in cost-to-go functions that are convex with respect to some state
variables (e.g., the reservoir levels) and concave with respect to other state
variables (e.g., the spot price in the current stage).

To overcome this problem, the approach in the literature has been to discretize
the price process in order to model it using a Markovian policy graph like those
discussed in [Basic IV: Markov uncertainty](@ref).

However, recent work offers a way to include stagewise-dependent objective
uncertainty into the objective function of SDDP subproblems. Readers are
directed to the following works for an introduction:

 - Downward, A., Dowson, O., and Baucke, R. (2017). Stochastic dual dynamic
   programming with stagewise dependent objective uncertainty. Optimization
   Online. [link](http://www.optimization-online.org/DB_HTML/2018/02/6454.html)

 - Dowson, O. PhD Thesis. University of Auckland, 2018. [link](https://researchspace.auckland.ac.nz/handle/2292/37700)

The method discussed in the above works introduces the concept of an _objective
state_ into SDDP. Unlike normal state variables in SDDP (e.g., the volume of
water in the reservoir), the cost-to-go function is _concave_ with respect to
the objective states. Thus, the method builds an outer approximation of the
cost-to-go function in the normal state-space, and an inner approximation of the
cost-to-go function in the objective state-space.

!!! warn
    Support for objective states in `SDDP.jl` is experimental. Models are
    considerably more computational intensive, the interface is less
    user-friendly, and there are [subtle gotchas to be aware of](@ref
    objective_state_warnings). Only use this if you have read and understood the
    theory behind the method.

## One-dimensional objective states

Let's assume that the fuel cost is not fixed, but instead evolves according to a
multiplicative auto-regressive process: `fuel_cost[t] = ω * fuel_cost[t-1]`,
where `ω` is drawn from the sample space `[0.75, 0.9, 1.1, 1.25]` with equal
probability.

An objective state can be added to a subproblem using the
[`Kokako.add_objective_state`](@ref) function. This can only be called once per
subproblem. If you want to add a multi-dimensional objective state, read
[Multi-dimensional objective states](@ref). [`Kokako.add_objective_state`](@ref)
takes a number of keyword arguments. The two required ones are

 - `initial_value`: the value of the objective state at the root node of the
   policy graph (i.e., identical to the `initial_value` when defining normal
   state variables.

 - `lipschitz`: the Lipschitz constant of the cost-to-go function with respect
   to the objective state. In other words, this value is the maximum change in
   the cost-to-go function _at any point in the state space_, given a one-unit
   change in the objective state.

There are also two optional keyword arguments: `lower_bound` and `upper_bound`,
which give SDDP.jl hints (importantly, not constraints) about the domain of the
objective state. Setting these bounds appropriately can improve the speed of
convergence.

Finally, [`Kokako.add_objective_state`](@ref) requires an update function. This
function takes two arguments. The first is the incoming value of the objective
state, and the second is the realization of the stagewise-independent noise term
(set using [`Kokako.parameterize`](@ref)). The function should return the value
of the objective state to be used in the current subproblem.

This connection with the stagewise-independent noise term means that
[`Kokako.parameterize`](@ref) _must_ be called in a subproblem that defines an
objective state. Inside [`Kokako.parameterize`](@ref), the value of the objective
state to be used in the current subproblem (i.e., after the update function),
can be queried using [`Kokako.objective_state`](@ref).

Here is the full model with the objective state.

```jldoctest intermediate_price
using Kokako, GLPK

model = Kokako.LinearPolicyGraph(
            stages = 3, sense = :Min, lower_bound = 0.0,
            optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, t
    @variable(subproblem, 0 <= volume <= 200, Kokako.State, initial_value = 200)
    @variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        inflow
    end)
    @constraints(subproblem, begin
        volume.out == volume.in + inflow - hydro_generation - hydro_spill
        demand_constraint, thermal_generation + hydro_generation == 150.0
    end)

    ###
    ### Add an objective state. ω will be the same value that is called in
    ### `Kokako.parameterize`.
    ###

    Kokako.add_objective_state(
            subproblem, initial_value = 50.0, lipschitz = 10_000.0,
            lower_bound = 50.0, upper_bound = 150.0) do fuel_cost, ω
        return ω.fuel * fuel_cost
    end

    ###
    ### Create the cartesian product of a multi-dimensional random variable.
    ###

    Ω = [
        (fuel = f, inflow = w)
        for f in [0.75, 0.9, 1.1, 1.25] for w in [0.0, 50.0, 100.0]
    ]

    Kokako.parameterize(subproblem, Ω) do ω
        ###
        ### Query the current fuel cost.
        ###

        fuel_cost = Kokako.objective_state(subproblem)

        @stageobjective(subproblem, fuel_cost * thermal_generation)
        JuMP.fix(inflow, ω.inflow)
    end
end

# output

A policy graph with 3 nodes.
 Node indices: 1, 2, 3
```

After creating our model, we can train and simulate as usual.

```jldoctest intermediate_price
Kokako.train(model, iteration_limit = 10, perform_numerical_stability_check=false)

simulations = Kokako.simulate(model, 1)

print("Finished training and simulating.")

# output

----------------------------------------------------
         SDDP.jl (c) Oscar Dowson, 2017-19

 Iteration   Simulation      Bound        Time (s)
         1   6.44531e+03   3.33333e+03   3.73000e-01
         2   7.56250e+03   4.53474e+03   3.76000e-01
         3   3.09375e+03   4.79404e+03   3.80000e-01
         4   2.06250e+03   4.79404e+03   3.83000e-01
         5   3.62399e+03   4.83141e+03   3.86000e-01
         6   7.24687e+03   4.86192e+03   3.90000e-01
         7   5.46875e+03   4.89025e+03   3.93000e-01
         8   1.68750e+03   4.89496e+03   4.11000e-01
         9   1.00781e+04   4.89496e+03   4.15000e-01
        10   1.05469e+04   4.89496e+03   4.19000e-01

Terminating training with status: iteration_limit
----------------------------------------------------
Finished training and simulating.
```

To demonstrate how the objective states are updated, consider the sequence of
noise observations:
```jldoctest intermediate_price; filter=r"\(fuel \= .+?\.0\)"
julia> [stage[:noise_term] for stage in simulations[1]]
3-element Array{NamedTuple{(:fuel, :inflow),Tuple{Float64,Float64}},1}:
 (fuel = 0.75, inflow = 0.0)
 (fuel = 0.9, inflow = 50.0)
 (fuel = 1.25, inflow = 50.0)
```

This, the fuel cost in the first stage should be `0.75 * 50 = 37.5`. The fuel
cost in the second stage should be `0.9 * 37.5 = 33.75`. The fuel cost in the
third stage should be `1.25 * 33.75 = 42.1875`.

To confirm this, the values of the objective state in a simulation can be
queried using the `:objective_state` key.

```jldoctest intermediate_price; filter=r"\d+?\.\d+"
julia> [stage[:objective_state] for stage in simulations[1]]
3-element Array{Float64,1}:
 37.5
 33.75
 42.1875
```

## Multi-dimensional objective states

You can construct multi-dimensional price processes using `NTuple`s. Just
replace every scalar value associated with the objective state by a tuple. For
example, `initial_value = 1.0` becomes `initial_value = (1.0, 2.0)`.

Here is an example:

```jldoctest intermediate_price
model = Kokako.LinearPolicyGraph(
            stages = 3, sense = :Min, lower_bound = 0.0,
            optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, t
    @variable(subproblem, 0 <= volume <= 200, Kokako.State, initial_value = 200)
    @variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        inflow
    end)
    @constraints(subproblem, begin
        volume.out == volume.in + inflow - hydro_generation - hydro_spill
        demand_constraint, thermal_generation + hydro_generation == 150.0
    end)

    Kokako.add_objective_state(
            subproblem, initial_value = (50.0, 50.0),
            lipschitz = (10_000.0, 10_000.0), lower_bound = (50.0, 50.0),
            upper_bound = (150.0, 150.0)) do fuel_cost, ω
        fuel_cost′ = fuel_cost[1] + 0.5 * (fuel_cost[1] - fuel_cost[2]) + ω.fuel
        return (fuel_cost′, fuel_cost[1])
    end

    Ω = [
        (fuel = f, inflow = w)
        for f in [-10.0, -5.0, 5.0, 10.0] for w in [0.0, 50.0, 100.0]
    ]

    Kokako.parameterize(subproblem, Ω) do ω
        (fuel_cost, fuel_cost_old) = Kokako.objective_state(subproblem)
        @stageobjective(subproblem, fuel_cost * thermal_generation)
        JuMP.fix(inflow, ω.inflow)
    end
end

Kokako.train(model, iteration_limit = 10, perform_numerical_stability_check=false)

simulations = Kokako.simulate(model, 1)

print("Finished training and simulating.")

# output

----------------------------------------------------
         SDDP.jl (c) Oscar Dowson, 2017-19

 Iteration   Simulation      Bound        Time (s)
         1   7.87500e+03   2.88826e+03   3.55000e-01
         2   4.25000e+03   4.08745e+03   3.59000e-01
         3   6.00000e+03   4.08745e+03   3.62000e-01
         4   1.62500e+03   4.08745e+03   3.66000e-01
         5   7.00000e+03   4.29941e+03   3.70000e-01
         6   8.12500e+02   4.36703e+03   3.73000e-01
         7   1.56250e+03   4.36703e+03   3.77000e-01
         8   7.37500e+03   4.78689e+03   3.80000e-01
         9   7.50000e+02   4.78689e+03   3.85000e-01
        10   5.39216e+03   4.84133e+03   3.88000e-01

Terminating training with status: iteration_limit
----------------------------------------------------
Finished training and simulating.
```

This time, since our objective state is two-dimensional, the objective states
are tuples with two elements:

```jldoctest intermediate_price; filter=r"\(.+?\)"
julia> [stage[:objective_state] for stage in simulations[1]]
3-element Array{Tuple{Float64,Float64},1}:
 (55.0, 50.0)
 (52.5, 55.0)
 (56.25, 52.5)
```

## [Warnings](@id objective_state_warnings)

There are number of things to be aware of when using objective states.

 - The key assumption is that price is independent of the states and actions in
    the model.

    That means that the price cannot appear in any `@constraint`s. Nor can you
    use any `@variable`s in the update function.

 - Choosing an appropriate Lipschitz constant is difficult.

    The points discussed in [Choosing an initial bound](@ref) are relevant. The
    Lipschitz constant should not be chosen as large as possible (since this
    will help with convergence and the numerical issues discussed above), but if
    chosen to small, it may cut of the feasible region and lead to a sub-optimal
    solution.

 - You need to ensure that the cost-to-go function is concave with respect to
    the objective state _before_ the update.

    $V(x, y) = \min\{Y(y)' x \;|\; Ax \ge b\}$

    If the update function is linear, this is always the case. In some
    situations, the update function can be nonlinear (e.g., multiplicative as we
    have above). In general, placing constraints on the price (e.g.,
    `clamp(price, 0, 1)`) will destroy concavity. [Caveat
    emptor](https://en.wikipedia.org/wiki/Caveat_emptor). It's up to you if this
    is a problem. If it isn't you'll get a good heuristic with no guarantee of
    global optimality.

In the next tutorial, [Intermediate V: performance](@ref), we discuss how to
improve the computational performance of `SDDP.jl` models.
