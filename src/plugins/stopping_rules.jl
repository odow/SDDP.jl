"""
    AbstractStoppingRule

The abstract type for the stopping-rule interface.

You need to define the following methods:
 - Kokako.stopping_rule_status
 - Kokako.convergence_test
"""
abstract type AbstractStoppingRule end

function stopping_rule_status(::AbstractStoppingRule)
    error("Not yet implemented.")
end

function convergence_test(graph::PolicyGraph, ::AbstractStoppingRule)
    error("Not yet implemented.")
end

function convergence_test(graph::PolicyGraph,
                          stopping_rules::Vector{AbstractStoppingRule})
    for stopping_rule in stopping_rules
        if convergence_test(graph, stopping_rule)
            return true, stopping_rule_status(stopping_rule)
        end
    end
    return false, :not_solved
end

# ========================= Statistical Stopping Rule ======================== #

struct Statistical <: AbstractStoppingRule end
stopping_rule_status(::Statistical) = :statistical
convergence_test(graph, ::Statistical) = error("Not yet implemented.")

# ======================= Bound-stalling Stopping Rule ======================= #

struct BoundStalling <: AbstractStoppingRule end
stopping_rule_status(::BoundStalling) = :bound_stalling
convergence_test(graph, ::BoundStalling) = error("Not yet implemented.")
