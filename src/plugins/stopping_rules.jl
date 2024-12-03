#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.
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
        num_replications::Int,
        iteration_period::Int = 1,
        z_score::Float64 = 1.96,
        verbose::Bool = true,
        disable_warning::Bool = false,
    )

Perform an in-sample Monte Carlo simulation of the policy with
`num_replications` replications every `iteration_period`s and terminate if the
deterministic bound (lower if minimizing) falls into the confidence interval for
the mean of the simulated cost.

If `verbose = true`, print the confidence interval.

If `disable_warning = true`, disable the warning telling you not to use this
stopping rule (see below).

## Why this stopping rule is not good

This stopping rule is one of the most common stopping rules seen in the
literature. Don't follow the crowd. It is a poor choice for your model, and
should be rarely used. Instead, you should use the default stopping rule, or use
a fixed limit like a time or iteration limit.

To understand why this stopping rule is a bad idea, assume we have conducted
`num_replications` simulations and the objectives are in a vector
`objectives::Vector{Float64}`.

Our mean is `μ = mean(objectives)` and the half-width of the confidence interval
is `w = z_score * std(objectives) / sqrt(num_replications)`.

Many papers suggest terminating the algorithm once the deterministic bound
(lower if minimizing, upper if maximizing) is contained within the confidence
interval. That is, if `μ - w <= bound <= μ + w`. Even worse, some papers define
an optimization gap of `(μ + w) / bound` (if minimizing) or `(μ - w) / bound`
(if maximizing), and they terminate once the gap is less than a value like 1%.

Both of these approaches are misleading, and more often than not, they will
result in terminating with a sub-optimal policy that performs worse than
expected. There are two main reasons for this:

 1. The half-width depends on the number of replications. To reduce the
    computational cost, users are often tempted to choose a small number of
    replications. This increases the half-width and makes it more likely that
    the algorithm will stop early. But if we choose a large number of
    replications, then the computational cost is high, and we would have been
    better off to run a fixed number of iterations and use that computational
    time to run extra training iterations.
 2. The confidence interval assumes that the simulated values are normally
    distributed. In infinite horizon models, this is almost never the case. The
    distribution is usually closer to exponential or log-normal.

There is a third, more technical reason which relates to the conditional
dependence of constructing multiple confidence intervals.

The default value of `z_score = 1.96` corresponds to a 95% confidence
interval. You should interpret the interval as "if we re-run this simulation
100 times, then the true mean will lie in the confidence interval 95 times
out of 100." But if the bound is within the confidence interval, then we know
the true mean cannot be better than the bound. Therfore, there is a more than
95% chance that the mean is within the interval.

A separate problem arises if we simulate, find that the bound is outside the
confidence interval, keep training, and then re-simulate to compute a new
confidence interval. Because we will terminate when the bound enters the
confidence interval, the repeated construction of a confidence interval means
that the unconditional probability that we terminate with a false positive is
larger than 5% (there are now more chances that the sample mean is optimistic
and that the confidence interval includes the bound but not the true mean). One
fix is to simulate with a sequentially increasing number of replicates, so that
the unconditional probability stays at 95%, but this runs into the problem of
computational cost. For more information on sequential sampling, see, for
example, Güzin Bayraksan, David P. Morton, (2011) A Sequential Sampling
Procedure for Stochastic Programming. Operations Research 59(4):898-913.
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
        disable_warning::Bool = false,
    )
        if !disable_warning
            @warn(
                "Are you really sure you want to use this stopping rule? " *
                "Read why we don't recommend it by typing `? " *
                "SDDP.Statistical` into the REPL to read the docstring.\n\n" *
                "If you understand what you are doing, you can disable this " *
                "warning with `SDDP.Statistical(; disable_warning = true)`",
            )
        end
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
    else
        @assert graph.objective_sense == MOI.MAX_SENSE
        return current_bound <= sample_mean + sample_ci
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
    distance_tol::Float64
    bound_tol::Float64
end

"""
    SimulationStoppingRule(;
        sampling_scheme::AbstractSamplingScheme = SDDP.InSampleMonteCarlo(),
        replications::Int = -1,
        period::Int = -1,
        distance_tol::Float64 = 1e-2,
        bound_tol::Float64 = 1e-4,
    )

Terminate the algorithm using a mix of heuristics. Unless you know otherwise,
this is typically a good default.

## Termination criteria

First, we check that the deterministic bound has stabilized. That is, over the
last five iterations, the deterministic bound has changed by less than an
absolute or relative tolerance of `bound_tol`.

Then, if we have not done one in the last `period` iterations, we perform a
primal simulation of the policy using `replications` out-of-sample realizations
from `sampling_scheme`. The realizations are stored and re-used in each
simulation. From each simulation, we record the value of the stage objective.
We terminate the policy if each of the trajectories in two consecutive
simulations differ by less than `distance_tol`.

By default, `replications` and `period` are `-1`, and SDDP.jl will guess good
values for these. Over-ride the default behavior by setting an appropriate
value.

## Example

```julia
SDDP.train(model; stopping_rules = [SimulationStoppingRule()])
```
"""
function SimulationStoppingRule(;
    sampling_scheme::AbstractSamplingScheme = InSampleMonteCarlo(),
    replications::Int = -1,
    period::Int = -1,
    distance_tol::Float64 = 1e-2,
    bound_tol::Float64 = 1e-4,
)
    cached_sampling_scheme =
        PSRSamplingScheme(replications; sampling_scheme = sampling_scheme)
    function simulator(model, N)
        cached_sampling_scheme.N = max(N, cached_sampling_scheme.N)
        scenarios = simulate(
            model,
            N;
            sampling_scheme = cached_sampling_scheme,
            # TODO(odow): if we use Threaded() here, then we get out of order
            # between the simulations and the PSRSamplingScheme: it's not the
            # case that simulation `i` accesses sample `i`, because they might
            # happen out-of-order.
            parallel_scheme = Serial(),
        )
        # !!! info
        #     At one point, I tried adding the primal value of the state
        #     variables. But it didn't work for some models because of
        #     degeneracy, that is, the value of a state variable will oscillate
        #     between two equally optimal outcomes in subsequent iterations.
        #     So for now, I just use the stage objective and the bellman term.
        keys = [:stage_objective, :bellman_term]
        return map(scenarios) do scenario
            return [getindex.(scenario, k) for k in keys]
        end
    end
    return SimulationStoppingRule(
        simulator,
        replications,
        period,
        Any[],
        0,
        distance_tol,
        bound_tol,
    )
end

stopping_rule_status(::SimulationStoppingRule) = :simulation_stopping

function _compute_distance(x::Real, y::Real)
    if x ≈ y
        return 0.0
    end
    return abs(x - y) / max(1.0, abs(x), abs(y))
end

function _compute_distance(new_data::Vector, old_data::Vector)
    d = sum(_compute_distance(x, y)^2 for (x, y) in zip(new_data, old_data))
    return sqrt(d)
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
    if !isapprox(
        log[end].bound,
        log[end-5].bound;
        atol = rule.bound_tol,
        rtol = rule.bound_tol,
    )
        return false  # If the lower bound haven't stalled, keep going.
    end
    if length(log) - rule.last_iteration < _period(rule.period, length(log))
        return false  # Do at least rule.period iterations since the last trial
    end
    new_data = rule.simulator(model, rule.replications)
    distance = _compute_distance(new_data, rule.data)
    rule.data = new_data
    rule.last_iteration = length(log)
    return distance < rule.distance_tol
end

# ========================== FirstStageStoppingRule ========================== #

mutable struct FirstStageStoppingRule <: AbstractStoppingRule
    data::Vector{Any}
    atol::Float64
    iterations::Int
end

"""
    FirstStageStoppingRule(; atol::Float64 = 1e-3, iterations::Int = 50)

Terminate the algorithm when the outgoing values of the first-stage state
variables have not changed by more than `atol` for `iterations` number of
consecutive iterations.

## Example

```julia
SDDP.train(model; stopping_rules = [FirstStageStoppingRule()])
```
"""
function FirstStageStoppingRule(; atol::Float64 = 1e-3, iterations::Int = 50)
    return FirstStageStoppingRule(Any[], atol, iterations)
end

stopping_rule_status(::FirstStageStoppingRule) = :first_stage_stopping

function convergence_test(
    model::PolicyGraph{T},
    log::Vector{Log},
    rule::FirstStageStoppingRule,
) where {T}
    if length(model.root_children) != 1
        error(
            "FirstStageStoppingRule cannot be applied because first-stage is " *
            "not deterministic",
        )
    end
    node = model[model.root_children[1].term]
    if length(node.noise_terms) > 1
        error(
            "FirstStageStoppingRule cannot be applied because first-stage is " *
            "not deterministic",
        )
    end
    lock(node.lock)
    try
        set_incoming_state(node, model.initial_root_state)
        parameterize(node, first(node.noise_terms).term)
        optimize!(node.subproblem)
        state = get_outgoing_state(node)
        push!(rule.data, state)
        if length(rule.data) < rule.iterations
            return false
        end
        for i in 1:(rule.iterations-1), (k, v) in state
            if !isapprox(rule.data[end-i][k], v; atol = rule.atol)
                return false
            end
        end
        return true
    finally
        unlock(node.lock)
    end
end
