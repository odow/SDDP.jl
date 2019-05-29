#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

# ================================ risk_measures ============================= #

"""
    AbstractRiskMeasure

The abstract type for the risk measure interface.

You need to define the following methods:
 - [`SDDP.adjust_probability`](@ref)
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
 - [`SDDP.sample_scenario`](@ref)
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
function sample_scenario end

# ============================== bellman_functions =========================== #

"""
    AbstractBellmanFunction

The abstract type for the Bellman function interface.

You need to define the following methods:
 - [`SDDP.initialize_bellman_function`](@ref)
 - [`SDDP.refine_bellman_function`](@ref)
 - [`SDDP.bellman_term`](@ref)
"""
abstract type AbstractBellmanFunction end

"""
    initialize_bellman_function(
        ::Type{F}, graph::PolicyGraph{T}, node::Node{T}
        ) where {F<:AbstractBellmanFunction, T}

Return an instance of the Bellman function F for `node` in the policy graph
`graph`.
"""
function initialize_bellman_function end

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
function refine_bellman_function end

"""
    bellman_term(::AbstractBellmanFunction)

Return a JuMP expression representing the Bellman function.
"""
function bellman_term end

# =============================== stopping_rules ============================= #

"""
    AbstractStoppingRule

The abstract type for the stopping-rule interface.

You need to define the following methods:
 - [`SDDP.stopping_rule_status`](@ref)
 - [`SDDP.convergence_test`](@ref)
"""
abstract type AbstractStoppingRule end

"""
    stopping_rule_status(::AbstractStoppingRule)::Symbol

Return a symbol describing the stopping rule.
"""
function stopping_rule_status end

"""
    convergence_test(model::PolicyGraph, log::Vector{Log}, ::AbstractStoppingRule)::Bool

Return a `Bool` indicating if the algorithm should terminate the training.
"""
function convergence_test end

function convergence_test(graph::PolicyGraph, log::Vector{Log},
                          stopping_rules::Vector{AbstractStoppingRule})
    for stopping_rule in stopping_rules
        if convergence_test(graph, log, stopping_rule)
            return true, stopping_rule_status(stopping_rule)
        end
    end
    return false, :not_solved
end
