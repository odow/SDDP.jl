#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

# ================================ risk_measures ============================= #

"""
    AbstractRiskMeasure

The abstract type for the risk measure interface.

You need to define the following methods:
 - Kokako.adjust_probability
"""
abstract type AbstractRiskMeasure end

"""
    adjust_probability(measure::Expectation
                       risk_adjusted_probability::Vector{Float64},
                       original_probability::Vector{Float64},
                       noise_support::Vector{Noise{T}},
                       objective_realizations::Vector{Float64},
                       is_minimization::Bool) where T
"""
function adjust_probability end

# ============================== sampling_schemes ============================ #

"""
    AbstractSamplingScheme

The abstract type for the sampling-scheme interface.

You need to define the following methods:
 - Kokako.sample_scenario
"""
abstract type AbstractSamplingScheme end

"""
    sample_scenario(graph::PolicyGraph{T}, ::AbstractSamplingScheme) where T

Sample a scenario from the policy graph `graph` based on the sampling scheme.

Returns `::Tuple{Vector{Tuple{T, <:Any}}, Bool}`, where the first element is the
scenario, and the second element is a Boolean flag indicating if the scenario
was terminated due to the detection of a cycle.

The scenario is a list of tuples (type `Vector{Tuple{T, <:Any}}`) where the
first component of each tuple is the index of the node, and the second component
is the stagewise-independent noise term observed in that node.
"""
function sample_scenario(graph::PolicyGraph{T},
                         sampling_scheme::AbstractSamplingScheme) where T
    error("You need to overload the function Kokako.sample_scenario for the " *
          "sampling scheme (sampling_scheme).")
end

# ============================== bellman_functions =========================== #

"""
    AbstractBellmanFunction

The abstract type for the Bellman function interface.

You need to define the following methods:
 - Kokako.initialize_bellman_function
 - Kokako.refine_bellman_function
 - Kokako.bellman_term
"""
abstract type AbstractBellmanFunction end

"""
    initialize_bellman_function(::Type{F}, graph::PolicyGraph{T}, node::Node{T}
                                    ) where {F<:AbstractBellmanFunction, T}

Return an instance of the Bellman function F for `node` in the policy graph
`graph`.
"""
function initialize_bellman_function(
        ::Type{F}, graph::PolicyGraph{T}, node::Node{T}
            ) where {F<:AbstractBellmanFunction, T}
    error("Overload the function Kokako.initialize_bellman_function for $(F).")
end

"""
    refine_bellman_function(graph::PolicyGraph{T},
                            node::Node{T},
                            bellman_function::AbstractBellmanFunction,
                            risk_measure::AbstractRiskMeasure,
                            state::Dict{Symbol, Float64},
                            dual_variables::Vector{Dict{Symbol, Float64}},
                            noise_supports::Vector{<:Noise},
                            original_probability::Vector{Float64},
                            objective_realizations::Vector{Float64}
                                ) where T
"""
function refine_bellman_function(graph::PolicyGraph{T},
                                 node::Node{T},
                                 bellman_function::AbstractBellmanFunction,
                                 risk_measure::AbstractRiskMeasure,
                                 outgoing_state::Dict{Symbol, Float64},
                                 dual_variables::Vector{Dict{Symbol, Float64}},
                                 noise_supports::Vector,
                                 original_probability::Vector{Float64},
                                 objective_realizations::Vector{Float64}
                                     ) where T
    error("Kokako.refine_bellman_function not implemented for " *
          "$(bellman_function).")
end

"""
    bellman_term(::AbstractBellmanFunction)

Return a JuMP expression representing the Bellman function.
"""
function bellman_term(bellman::AbstractBellmanFunction)
    error("Kokako.bellman term not implemented for $(bellman).")
end

# =============================== stopping_rules ============================= #

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
