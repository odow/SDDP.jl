"""
    AbstractStoppingRule

The abstract type for the stopping-rule interface.

You need to define the following methods:
 - Kokako.stopping_rule_status
 - Kokako.convergence_test
"""
abstract type AbstractStoppingRule end

struct Log
    iteration::Int
    bound::Float64
    simulation_value::Float64
    time::Float64
end

function stopping_rule_status(stopping_rule::AbstractStoppingRule)
    error("You need to overload the function Kokako.stopping_rule_status for " *
          "the stopping rule (stopping_rule).")
end

function convergence_test(
    graph::PolicyGraph, log::Vector{Log}, stopping_rule::AbstractStoppingRule)
    error("You need to overload the function Kokako.convergence_test for the " *
          "stopping rule (stopping_rule).")
end

function convergence_test(graph::PolicyGraph,
                          log::Vector{Log},
                          stopping_rules::Vector{AbstractStoppingRule})
    for stopping_rule in stopping_rules
        if convergence_test(graph, log, stopping_rule)
            return true, stopping_rule_status(stopping_rule)
        end
    end
    return false, :not_solved
end

# ======================= Iteration Limit Stopping Rule ====================== #

mutable struct IterationLimit <: AbstractStoppingRule
    limit::Int
end

stopping_rule_status(::IterationLimit) = :iteration_limit

function convergence_test(
        graph::PolicyGraph, log::Vector{Log}, rule::IterationLimit)
    return log[end].iteration >= rule.limit
end

# ========================= Time Limit Stopping Rule ========================= #

mutable struct TimeLimit <: AbstractStoppingRule
    limit::Float64
end

stopping_rule_status(::TimeLimit) = :time_limit

function convergence_test(graph::PolicyGraph, log::Vector{Log}, rule::TimeLimit)
    return log[end].time >= rule.limit
end

# ========================= Statistical Stopping Rule ======================== #

# struct Statistical <: AbstractStoppingRule
    # number_replications::Int
# end
#
# stopping_rule_status(::Statistical) = :statistical
# function convergence_test(graph, ::Statistical)
#
# end

# ======================= Bound-stalling Stopping Rule ======================= #

struct BoundStalling <: AbstractStoppingRule
    num_previous_iterations::Int
    tolerance::Float64
end

stopping_rule_status(::BoundStalling) = :bound_stalling

function convergence_test(graph::PolicyGraph{T}, log::Vector{Log},
                          rule::BoundStalling) where T
    if length(log) < rule.num_previous_iterations + 1
        return false
    end
    for iteration in 1:rule.num_previous_iterations
        idx = length(log) - iteration
        if abs(log[idx].bound - log[idx + 1].bound) > rule.tolerance
            return false
        end
    end
    return true
end
