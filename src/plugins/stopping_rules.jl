#  Copyright 2017-20, Oscar Dowson.
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

function convergence_test(graph::PolicyGraph, log::Vector{Log}, rule::IterationLimit)
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

function convergence_test(graph::PolicyGraph, log::Vector{Log}, rule::TimeLimit)
    return log[end].time >= rule.limit
end

# ========================= Statistical Stopping Rule ======================== #

"""
    Statistical(; num_replications, iteration_period = 1, z_score = 1.96,
                verbose = true)

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

function convergence_test(graph::PolicyGraph, log::Vector{Log}, rule::Statistical)
    if length(log) % rule.iteration_period != 0
        # Only run this convergence test every rule.iteration_period iterations.
        return false
    end
    results = simulate(graph, rule.num_replications)
    objectives = map(simulation -> sum(s[:stage_objective] for s in simulation), results)
    sample_mean = Statistics.mean(objectives)
    sample_ci = rule.z_score * Statistics.std(objectives) / sqrt(rule.num_replications)
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
        #If sense is none of the above for some awkward reason, return to previous criteria
        return sample_mean - sample_ci <= current_bound <= sample_mean + sample_ci
    end

end

# ======================= Bound-stalling Stopping Rule ======================= #

"""
    BoundStalling(num_previous_iterations::Int, tolerance::Float64)

Teriminate the algorithm once the deterministic bound (lower if minimizing,
upper if maximizing) fails to improve by more than `tolerance` in absolute terms
for more than `num_previous_iterations` consecutve iterations.
"""
struct BoundStalling <: AbstractStoppingRule
    num_previous_iterations::Int
    tolerance::Float64
end

stopping_rule_status(::BoundStalling) = :bound_stalling

function convergence_test(
    graph::PolicyGraph{T},
    log::Vector{Log},
    rule::BoundStalling,
) where {T}
    if length(log) < rule.num_previous_iterations + 1
        return false
    end
    for iteration = 1:rule.num_previous_iterations
        idx = length(log) - iteration
        if abs(log[idx].bound - log[idx+1].bound) > rule.tolerance
            return false
        end
    end
    return true
end
