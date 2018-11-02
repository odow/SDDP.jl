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

# ============================== SDDP.AverageCut ===============================

mutable struct Cut
    intercept::Float64
    coefficients::Dict{Symbol, Float64}
    index
end

struct AverageCut <: AbstractBellmanFunction
    variable::JuMP.VariableRef
    cut_improvement_tolerance::Float64
    cuts::Vector{Cut}
end

struct BellmanFactory{T}
    args
    kwargs
    BellmanFactory{T}(args...; kwargs...) where T = new{T}(args, kwargs)
end

"""
    AverageCut(; lower_bound = -Inf, upper_bound = Inf)

The AverageCut Bellman function. Provide a lower_bound if minimizing, or an
upper_bound if maximizing.
"""
# function AverageCut(;
#         lower_bound = -Inf,
#         upper_bound = Inf,
#         cut_improvement_tolerance::Float64 = 0.0
#         )
function AverageCut(; kwargs...)
    return BellmanFactory{AverageCut}(; kwargs...)
end

function initialize_bellman_function(factory::BellmanFactory{AverageCut},
                                     graph::PolicyGraph{T},
                                     node::Node{T}) where T
    lower_bound, upper_bound = -Inf, Inf
    cut_improvement_tolerance = 0.0
    if length(factory.args) > 0
        error("Positional arguments $(factory.args) ignored in AverageCut.")
    end
    for (kw, value) in factory.kwargs
        if kw == :lower_bound
            lower_bound = value
        elseif kw == :upper_bound
            upper_bound = value
        elseif kw == :cut_improvement_tolerance
            if value < 0
                error("Cut cut_improvement_tolerance must be > 0.")
            end
            cut_improvement_tolerance = value
        else
            error("Keyword $(kw) not recognised as argument to AverageCut.")
        end
    end
    bellman_variable = if length(node.children) > 0
        @variable(node.subproblem,
                  lower_bound = lower_bound, upper_bound = upper_bound)
    else
        @variable(node.subproblem, lower_bound = 0, upper_bound = 0)
    end
    return AverageCut(bellman_variable, cut_improvement_tolerance, Cut[])
end

bellman_term(bellman::AverageCut) = bellman.variable

function refine_bellman_function(graph::PolicyGraph{T},
                                 node::Node{T},
                                 bellman_function::AverageCut,
                                 risk_measure::AbstractRiskMeasure,
                                 outgoing_state::Dict{Symbol, Float64},
                                 dual_variables::Vector{Dict{Symbol, Float64}},
                                 noise_supports::Vector,
                                 original_probability::Vector{Float64},
                                 objective_realizations::Vector{Float64}
                                     ) where T
    is_minimization = graph.objective_sense == MOI.MinSense
    risk_adjusted_probability = similar(original_probability)
    adjust_probability(risk_measure,
                       risk_adjusted_probability,
                       original_probability,
                       noise_supports,
                       objective_realizations,
                       is_minimization)
    # Initialize average cut coefficients.
    intercept = 0.0
    coefficients = Dict{Symbol, Float64}()
    for state in keys(outgoing_state)
        coefficients[state] = 0.0
    end
    # Gather up coefficients for cut calculation.
    # β = F[λ]
    # α = F[θ] - βᵀ ̄x'
    # θ ≥ α + βᵀ x'
    for (objective, dual, prob) in zip(objective_realizations, dual_variables,
                                       risk_adjusted_probability)
        intercept += prob * objective
        for state in keys(outgoing_state)
            coefficients[state] += prob * dual[state]
        end
    end
    # Height of the cut at outgoing_state. We cache the value here for the
    # tolerance check that happens later.
    current_height = intercept
    # Calculate the intercept of the cut.
    for (name, value) in outgoing_state
        intercept -= coefficients[name] * value
    end

    # A structure to hold information about the cut. The third argument is
    # `nothing` because we haven't added it to the model.
    cut = Cut(intercept, coefficients, nothing)

    # Test whether we should add the new cut to the subproblem. We do this now
    # before collating the intercept to avoid twice the work.
    cut_is_an_improvement = if bellman_function.cut_improvement_tolerance > 0.0
        abs(JuMP.objective_value(node.subproblem) - current_height) >
            bellman_function.cut_improvement_tolerance
    else
        true
    end

    if cut_is_an_improvement
        index = if is_minimization
            @constraint(node.subproblem, bellman_function.variable >=
                intercept + sum(coefficients[name] * state.out
                    for (name, state) in node.states))
        else
            @constraint(node.subproblem, bellman_function.variable <=
                intercept + sum(coefficients[name] * state.out
                    for (name, state) in node.states))
        end
        # Store the index of the cut.
        cut.index = index
    end
    # Store the cut in the Bellman function.
    push!(bellman_function.cuts, cut)
    return
end
