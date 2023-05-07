#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

# ======================= Iteration Limit Stopping Rule ====================== #

"""
    IterationLimit(limit::Int)

Teriminate the algorithm after `limit` number of iterations.
"""
mutable struct IterationLimit <: AbstractStoppingRule
    limit::Int
end

stopping_rule_status(::IterationLimit) = :iteration_limit

function convergence_test(::PolicyGraph, log::Vector{Log}, rule::IterationLimit)
    return log[end].iteration >= rule.limit
end

# ========================= Time Limit Stopping Rule ========================= #

"""
    TimeLimit(limit::Float64)

Teriminate the algorithm after `limit` seconds of computation.
"""
mutable struct TimeLimit <: AbstractStoppingRule
    limit::Float64
end

stopping_rule_status(::TimeLimit) = :time_limit

function convergence_test(::PolicyGraph, log::Vector{Log}, rule::TimeLimit)
    return log[end].time >= rule.limit
end

# ========================= Statistical Stopping Rule ======================== #

"""
    Statistical(;
        num_replications,
        iteration_period = 1,
        z_score = 1.96,
        verbose = true,
    )

Perform an in-sample Monte Carlo simulation of the policy with
`num_replications` replications every `iteration_period`s. Terminate if the
deterministic bound (lower if minimizing) calls into the confidence interval for
the mean of the simulated cost. If `verbose = true`, print the confidence
interval.

Note that this tests assumes that the simulated values are normally distributed.
In infinite horizon models, this is almost never the case. The distribution is
usually closer to exponential or log-normal.
"""
struct Statistical <: AbstractStoppingRule
    num_replications::Int
    iteration_period::Int
    z_score::Float64
    verbose::Bool
    function Statistical(;
        num_replications,
        iteration_period = 1,
        z_score = 1.96,
        verbose = true,
    )
        return new(num_replications, iteration_period, z_score, verbose)
    end
end

stopping_rule_status(::Statistical) = :statistical

function convergence_test(
    graph::PolicyGraph,
    log::Vector{Log},
    rule::Statistical,
)
    if length(log) % rule.iteration_period != 0
        # Only run this convergence test every rule.iteration_period iterations.
        return false
    end
    results = simulate(graph, rule.num_replications)
    objectives =
        map(simulation -> sum(s[:stage_objective] for s in simulation), results)
    sample_mean = Statistics.mean(objectives)
    sample_ci =
        rule.z_score * Statistics.std(objectives) / sqrt(rule.num_replications)
    if rule.verbose
        println(
            "Simulated policy value: [",
            print_value(sample_mean - sample_ci),
            ", ",
            print_value(sample_mean + sample_ci),
            "]",
        )
    end
    current_bound = log[end].bound

    if graph.objective_sense == MOI.MIN_SENSE
        return sample_mean - sample_ci <= current_bound
    elseif graph.objective_sense == MOI.MAX_SENSE
        return current_bound <= sample_mean + sample_ci
    else
        # If sense is none of the above for some awkward reason, return to
        # previous criteria
        return sample_mean - sample_ci <=
               current_bound <=
               sample_mean + sample_ci
    end
end

# ======================= Bound-stalling Stopping Rule ======================= #

"""
    BoundStalling(num_previous_iterations::Int, tolerance::Float64)

Teriminate the algorithm once the deterministic bound (lower if minimizing,
upper if maximizing) fails to improve by more than `tolerance` in absolute terms
for more than `num_previous_iterations` consecutve iterations, provided it has
improved relative to the bound after the first iteration.

Checking for an improvement relative to the first iteration avoids early
termination in a situation where the bound fails to improve for the first `N`
iterations. This frequently happens in models with a large number of stages,
where it takes time for the cuts to propogate backward enough to modify the
bound of the root node.
"""
struct BoundStalling <: AbstractStoppingRule
    num_previous_iterations::Int
    tolerance::Float64
end

stopping_rule_status(::BoundStalling) = :bound_stalling

function convergence_test(
    ::PolicyGraph{T},
    log::Vector{Log},
    rule::BoundStalling,
) where {T}
    if length(log) < rule.num_previous_iterations + 1
        return false
    end
    # No change in the bound. There are three possibilities:
    #  1) we haven't added enough cuts
    #  2) the problem was deterministic or myopic
    #  3) there were existing cuts
    existing_solves = log[1].total_solves > log[end].total_solves / length(log)
    if !existing_solves && isapprox(log[1].bound, log[end].bound; atol = 1e-6)
        return all(l -> isapprox(l.bound, l.simulation_value; atol = 1e-6), log)
    end
    for i in 1:rule.num_previous_iterations
        if abs(log[end-i].bound - log[end-i+1].bound) > rule.tolerance
            return false
        end
    end
    return true
end

"""
    StoppingChain(rules::AbstractStoppingRule...)

Terminate once all of the `rules` are statified.

This stopping rule short-circuits, so subsequent rules are only tested if the
previous pass.

## Examples

A stopping rule that runs 100 iterations, then checks for the bound stalling:
```julia
StoppingChain(IterationLimit(100), BoundStalling(5, 0.1))
```
"""
struct StoppingChain <: AbstractStoppingRule
    rules::Vector{AbstractStoppingRule}

    function StoppingChain(rules::AbstractStoppingRule...)
        return new(collect(rules))
    end
end

function stopping_rule_status(rule::StoppingChain)
    return Symbol(join(stopping_rule_status.(rule.rules), " ∧ "))
end

function convergence_test(
    graph::PolicyGraph,
    log::Vector{Log},
    chain::StoppingChain,
)
    for rule in chain.rules
        if !convergence_test(graph, log, rule)
            return false
        end
    end
    return true
end

# ========================== SimulationStoppingRule ========================== #

mutable struct SimulationStoppingRule{F} <: AbstractStoppingRule
    simulator::F
    replications::Int
    period::Int
    data::Vector{Any}
    last_iteration::Int
end

function _get_state_variable_value(key)
    return sp -> JuMP.value(JuMP.variable_by_name(sp, "$(key)_out"))
end

"""
    SimulationStoppingRule(;
        sampling_scheme::AbstractSamplingScheme = SDDP.InSampleMonteCarlo(),
        replications::Int = -1,
        period::Int = -1,
    )

Terminate the algorithm using a mix of heuristics. Unless you know otherwise,
this is typically a good default.

## Termination criteria

First, we check that the deterministic bound has stabilized. That is, over the
last five iterations, the deterministic bound has changed by less than an
absolute or relative tolerance of `1e-4`.

Then, if we have not done one in the last `period` iterations, we perform a
primal simulation of the policy using `replications` out-of-sample realizations
from `sampling_scheme`. The realizations are stored and re-used in each
simulation. From each simulation, we record the values of the stage objective
and state variables. We terminate the policy if each of the trajectories in two
consecutive simulations differ by less than 1%.

By default, `replications` and `period` are `-1`, and SDDP.jl will guess good
values for these. Over-ride the default behavior by setting an appropriate
value.

## Example

```julia
SDDP.train(model; sampling_schemes = [SimulationStoppingRule()])
```
"""
function SimulationStoppingRule(;
    sampling_scheme::AbstractSamplingScheme = SDDP.InSampleMonteCarlo(),
    replications::Int = -1,
    period::Int = -1,
)
    cached_sampling_scheme =
        SDDP.PSRSamplingScheme(replications; sampling_scheme = sampling_scheme)
    function simulator(model, N)
        cached_sampling_scheme.N = max(N, cached_sampling_scheme.N)
        states = collect(keys(model.initial_root_state))
        custom_recorders = Dict{Symbol,Function}(
            k => _get_state_variable_value(k) for k in states
        )
        scenarios = SDDP.simulate(
            model,
            N;
            sampling_scheme = cached_sampling_scheme,
            custom_recorders = custom_recorders,
        )
        return map(scenarios) do scenario
            return [getindex.(scenario, k) for k in vcat(states, :bellman_term)]
        end
    end
    return SimulationStoppingRule(simulator, replications, period, Any[], 0)
end

stopping_rule_status(::SimulationStoppingRule) = :simulation_stopping

function _compute_distance_inner(x::Float64, y::Float64)
    return abs(x - y) / max(1.0, abs(x), abs(y))
end

function _compute_distance_inner(x::SDDP.State, y::SDDP.State)
    return sqrt(
        _compute_distance_inner(x.in, y.in)^2 +
        _compute_distance_inner(x.out, y.out)^2,
    )
end

function _compute_distance_inner(sim_1::AbstractArray, sim_2::AbstractArray)
    d = [_compute_distance_inner(x, y)^2 for (x, y) in zip(sim_1, sim_2)]
    return sqrt(sum(d))
end

function _compute_distance(new_data::Vector, old_data::Vector)
    d = [_compute_distance_inner(x, y)^2 for (x, y) in zip(new_data, old_data)]
    return sqrt(sum(d))
end

function _period(period, iterations)
    if period != -1
        return period
    elseif iterations <= 100
        return 20
    elseif iterations <= 1_000
        return 100
    else
        return 500
    end
end

function convergence_test(
    model::PolicyGraph{T},
    log::Vector{Log},
    rule::SimulationStoppingRule,
) where {T}
    # Setup parameters based on the model.
    if rule.replications == -1
        rule.replications = min(100, _unique_paths(model))
    end
    if isempty(rule.data)
        # On the first iteration, run a simulation and keep going.
        rule.data = rule.simulator(model, rule.replications)
        rule.last_iteration = 0
        return false
    end
    if length(log) <= 5
        return false  # Always do at least 5 iterations.
    end
    if !isapprox(log[end].bound, log[end-5].bound; atol = 1e-4, rtol = 1e-4)
        return false  # If the lower bound haven't stalled, keep going.
    end
    if length(log) - rule.last_iteration < _period(rule.period, length(log))
        return false  # Do at least rule.period iterations since the last trial
    end
    new_data = rule.simulator(model, rule.replications)
    distance = _compute_distance(new_data, rule.data)
    rule.data = new_data
    rule.last_iteration = length(log)
    return distance < 0.01
end
