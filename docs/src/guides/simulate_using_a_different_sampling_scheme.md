# Simulate using a different sampling scheme

```@meta
DocTestSetup = quote
    using SDDP, GLPK
end
```

By default, [`SDDP.simulate`](@ref) will simulate the policy using the
distributions of noise terms that were defined when the model was created. We
call these _in-sample_ simulations. However, in general the _in-sample_
distributions are an approximation of some underlying probabiltiy model which
we term the _true process_. Therefore, `SDDP.jl` makes it easy to simulate the
policy using different probability distrutions.

To demonsrate the different ways of simulating the policy, we're going to use
the model from the tutorial [Basic IV: Markov uncertainty](@ref).

```jldoctest sampling_schemes
using SDDP, GLPK

Ω = [
    (inflow = 0.0, fuel_multiplier = 1.5),
    (inflow = 50.0, fuel_multiplier = 1.0),
    (inflow = 100.0, fuel_multiplier = 0.75)
]

model = SDDP.MarkovianPolicyGraph(
    transition_matrices = Array{Float64, 2}[
        [ 1.0 ]',
        [ 0.75 0.25 ],
        [ 0.75 0.25 ; 0.25 0.75 ]
    ],
    sense = :Min,
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer
) do subproblem, node
    # Unpack the stage and Markov index.
    t, markov_state = node
    # Define the state variable.
    @variable(subproblem, 0 <= volume <= 200, SDDP.State, initial_value = 200)
    # Define the control variables.
    @variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        inflow
    end)
    # Define the constraints
    @constraints(subproblem, begin
        volume.out == volume.in + inflow - hydro_generation - hydro_spill
        thermal_generation + hydro_generation == 150.0
    end)
    # Note how we can use `markov_state` to dispatch an `if` statement.
    probability = if markov_state == 1  # wet climate state
        [1/6, 1/3, 1/2]
    else  # dry climate state
        [1/2, 1/3, 1/6]
    end

    fuel_cost = [50.0, 100.0, 150.0]
    SDDP.parameterize(subproblem, Ω, probability) do ω
        JuMP.fix(inflow, ω.inflow)
        @stageobjective(subproblem,
            ω.fuel_multiplier * fuel_cost[t] * thermal_generation)
    end
end

SDDP.train(model; iteration_limit = 10, print_level = 0);

model

# output

A policy graph with 5 nodes.
 Node indices: (1, 1), (2, 1), (2, 2), (3, 1), (3, 2)
```

## In-sample Monte Carlo simulation

To simulate the policy using the data defined when `model` was created, use
[`SDDP.InSampleMonteCarlo`](@ref).

```julia sampling_schemes
simulations = SDDP.simulate(
    model, 20, sampling_scheme = SDDP.InSampleMonteCarlo()
)

sort(unique(
    [node[:noise_term] for simulation in simulations for node in simulation]
))

# output

3-element Array{NamedTuple{(:inflow, :fuel_multiplier),Tuple{Float64,Float64}},1}:
 (inflow = 0.0, fuel_multiplier = 1.5)
 (inflow = 50.0, fuel_multiplier = 1.0)
 (inflow = 100.0, fuel_multiplier = 0.75)
```

## Out-of-sample Monte Carlo simulation

Instead of using the _in-sample_ data, we can perform an _out-of-sample_
simulation of the policy using the [`SDDP.OutOfSampleMonteCarlo`](@ref)
sampling scheme.

For each node, the [`SDDP.OutOfSampleMonteCarlo`](@ref) needs to define a new
distribution for the transition probabilities between nodes in the policy graph,
and a new distribution for the stagewise independent noise terms.

!!! note
    The support of the distribution for the stagewise independent noise terms
    does not have to be the same as the in-sample distributions.

```jldoctest sampling_schemes
sampling_scheme = SDDP.OutOfSampleMonteCarlo(model) do node
    stage, markov_state = node
    if stage == 0
        # Called from the root node. Transition to (1, 1) with probability 1.0.
        # Only return the list of children, _not_ a list of noise terms.
        return [SDDP.Noise((1, 1), 1.0)]
    elseif stage == 3
        # Called from the final node. Return an empty list for the children,
        # and a single, deterministic realization for the noise terms.
        children = SDDP.Noise[]
        noise_terms = [SDDP.Noise((inflow = 75.0, fuel_multiplier = 1.2), 1.0)]
        return children, noise_terms
    else
        # Called from a normal node. Return the in-sample distribution for the
        # noise terms, but modify the transition probabilities so that the
        # Markov switching probability is now 50%.
        probability = markov_state == 1 ? [1/6, 1/3, 1/2] : [1/2, 1/3, 1/6]
        noise_terms = [SDDP.Noise(ω, p) for (ω, p) in zip(Ω, probability)]
        children = [
            SDDP.Noise((stage + 1, 1), 0.5), SDDP.Noise((stage + 1, 2), 0.5)
        ]
        return children, noise_terms
    end
end

simulations = SDDP.simulate(model, 1, sampling_scheme = sampling_scheme)

simulations[1][3][:noise_term]

# output

(inflow = 75.0, fuel_multiplier = 1.2)
```

Alternatively, if you only want to modify the stagewise independent noise terms,
pass `use_insample_transition = true`.

```jldoctest sampling_schemes
sampling_scheme = SDDP.OutOfSampleMonteCarlo(
    model, use_insample_transition = true
) do node
stage, markov_state = node
    if stage == 3
        # Called from the final node. Return a single, deterministic
        # realization for the noise terms. Don't return the children because we
        # use the in-sample data.
        return [SDDP.Noise((inflow = 65.0, fuel_multiplier = 1.1), 1.0)]
    else
        # Called from a normal node. Return the in-sample distribution for the
        # noise terms. Don't return the children because we use the in-sample
        # data.
        probability = markov_state == 1 ? [1/6, 1/3, 1/2] : [1/2, 1/3, 1/6]
        return [SDDP.Noise(ω, p) for (ω, p) in zip(Ω, probability)]
    end
end

simulations = SDDP.simulate(model, 1, sampling_scheme = sampling_scheme)

simulations[1][3][:noise_term]

# output

(inflow = 65.0, fuel_multiplier = 1.1)
```

## Historical simulation

Instead of performing a Monte Carlo simulation like the previous tutorials, we
may want to simulate one particular sequence of noise realizations. This
_historical_ simulation can also be conducted by passing a
[`SDDP.Historical`](@ref) sampling scheme to the `sampling_scheme` keyword of
the [`SDDP.simulate`](@ref) function.

We can confirm that the historical sequence of nodes was visited by querying
the `:node_index` key of the simulation results.

```jldoctest sampling_schemes
simulations = SDDP.simulate(
    model,
    sampling_scheme = SDDP.Historical([
        ((1, 1), Ω[1]),
        ((2, 2), Ω[3]),
        ((3, 1), Ω[2])
    ])
)

[stage[:node_index] for stage in simulations[1]]

# output

3-element Array{Tuple{Int64,Int64},1}:
 (1, 1)
 (2, 2)
 (3, 1)
```
