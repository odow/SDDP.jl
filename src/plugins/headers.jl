#  Copyright 2017-21, Oscar Dowson.
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
    adjust_probability(
        measure::Expectation
        risk_adjusted_probability::Vector{Float64},
        original_probability::Vector{Float64},
        noise_support::Vector{Noise{T}},
        objective_realizations::Vector{Float64},
        is_minimization::Bool,
    ) where {T}
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
    sample_scenario(graph::PolicyGraph{T}, ::AbstractSamplingScheme) where {T}

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
        ::Type{F},
        graph::PolicyGraph{T},
        node::Node{T},
    ) where {F<:AbstractBellmanFunction,T}

Return an instance of the Bellman function F for `node` in the policy graph
`graph`.
"""
function initialize_bellman_function end

"""
    refine_bellman_function(
        graph::PolicyGraph{T},
        node::Node{T},
        bellman_function::AbstractBellmanFunction,
        risk_measure::AbstractRiskMeasure,
        state::Dict{Symbol, Float64},
        dual_variables::Vector{Dict{Symbol,Float64}},
        noise_supports::Vector{<:Noise},
        original_probability::Vector{Float64},
        objective_realizations::Vector{Float64},
    ) where {T}
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
    convergence_test(
        model::PolicyGraph,
        log::Vector{Log},
        ::AbstractStoppingRule,
    )::Bool

Return a `Bool` indicating if the algorithm should terminate the training.
"""
function convergence_test(
    graph::PolicyGraph,
    log::Vector{Log},
    stopping_rules::Vector{AbstractStoppingRule},
)
    for stopping_rule in stopping_rules
        if convergence_test(graph, log, stopping_rule)
            return true, stopping_rule_status(stopping_rule)
        end
    end
    return false, :not_solved
end

# ============================== backward_samplers =========================== #

"""
    AbstractBackwardSamplingScheme

The abstract type for backward sampling scheme interface.

You need to define the following methods:
 - [`SDDP.sample_backward_noise_terms`](@ref)
"""
abstract type AbstractBackwardSamplingScheme end

"""
    sample_backward_noise_terms(
        backward_sampling_scheme::AbstractBackwardSamplingScheme,
        node::Node{T},
    )::Vector{Noise}

Returns a `Vector{Noise}` of noises sampled from `node.noise_terms` using
`backward_sampling_scheme`
"""
function sample_backward_noise_terms end

# =========================== duality_handlers =========================== #

"""
    AbstractDualityHandler

The abstract type for the duality handler interface.
"""
abstract type AbstractDualityHandler end

"""
    get_dual_solution(
        node::Node,
        duality_handler::AbstractDualityHandler,
    )::Tuple{Float64,Dict{Symbol,Float64}}

Returns a `Float64` for the objective of the dual solution, and a
`Dict{Symbol,Float64}` where the keys are the names of the state variables and
the values are the dual variables associated with the fishing constraint at
`node`.
"""
function get_dual_solution end

"""
    prepare_backward_pass(
        node::Node,
        handler::AbstractDualityHandler,
        options::Options
    )

Performs any setup needed by the duality handler prior to the backward pass.

Returns a function that, when called with no arguments, undoes the setup.
"""
function prepare_backward_pass(::Node, ::AbstractDualityHandler, ::Options)
    return () -> nothing
end

# ============================= parallel schemes ============================= #

"""
    AbstractParallelScheme

Abstract type for different parallelism schemes.
"""
abstract type AbstractParallelScheme end

"""
    master_loop(
        ::AbstractParallelScheme,
        model::PolicyGraph{T},
        options::Options,
    )::Symbol where {T}

The solve loop of the SDDP algorithm. Returns a symbol corresponding to the
termination status.
"""
function master_loop end

"""
    _simulate(
        model::PolicyGraph,
        ::AbstractParallelScheme,
        number_replications::Int,
        variables::Vector{Symbol};
        kwargs...,
    )

Simulate the policy using the parallel scheme.
"""
function _simulate end

# ============================= forward pass ============================= #

"""
    AbstractForwardPass

Abstract type for different forward passes.
"""
abstract type AbstractForwardPass end

"""
    forward_pass(model::PolicyGraph, options::Options, ::AbstractForwardPass)

Return a forward pass as a named tuple with the following fields:

    (
        ;scenario_path,
        sampled_states,
        objective_states,
        belief_states,
        cumulative_value,
    )

See [`DefaultForwardPass`](@ref) for details.
"""
function forward_pass end
