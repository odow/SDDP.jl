#  Copyright 2017-21, Oscar Dowson.
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
