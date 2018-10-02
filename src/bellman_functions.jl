"""
    AbstractBellmanFunction

The abstract type for the Bellman function interface.

You need to define the following methods:
 - Kokako.initialize_bellman_function
 - Kokako.refine_bellman_function
 - Kokako.bellman_term
 - JSON.lower(bellman::AbstractBellmanFunction)
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
    error("Kokako.refine_bellman_function not implemented for $(bellman_function).")
end

"""
    bellman_term(::AbstractBellmanFunction)

Return a JuMP expression representing the Bellman function.
"""
function bellman_term(bellman::AbstractBellmanFunction)
    error("Kokako.bellman term not implemented for $(bellman).")
end

struct Cut
    intercept::Float64
    coefficients::Dict{Symbol, Float64}
end

function JSON.lower(cut::Cut)
    return Dict("intercept" => cut.intercept, "coefficients" => cut.coefficients)
end

"""
    write_bellman_to_file(policy_graph::PolicyGraph{T}, filename::String) where T

Save the Bellman function to `filename` in JSON format.
"""
function write_bellman_to_file(policy_graph::PolicyGraph{T},
                               filename::String) where T
    bellman = Dict{T, Any}()
    for (node_index, node) in policy_graph.nodes
        bellman[node_index] = Dict(
            "type" => typeof(node.bellman_function),
            "function" => node.bellman_function
        )
    end
    open(filename, "w") do io
        println(io, JSON.json(bellman))
    end
end

# ============================== SDDP.AverageCut ===============================

struct AverageCut <: AbstractBellmanFunction
    variable::JuMP.VariableRef
    cuts::Vector{Cut}
end

JSON.lower(bellman::AverageCut) = bellman.cuts

function initialize_bellman_function(
        ::Type{AverageCut}, graph::PolicyGraph{T}, node::Node{T}) where T
    bellman_variable = if length(node.children) > 0
        @variable(node.subproblem, lower_bound=-1000, upper_bound=1000)
    else
        @variable(node.subproblem, lower_bound=0, upper_bound=0)
    end
    return AverageCut(bellman_variable, Cut[])
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
    is_minimization = JuMP.objective_sense(node.subproblem) == :Min
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
        for (state, coefficient) in dual
            coefficients[state] += prob * coefficient
        end
    end
    for (name, value) in outgoing_state
        intercept -= coefficients[name] * value
    end
    # Add the cut to the subproblem.
    if is_minimization
        @constraint(node.subproblem, bellman_function.variable >=
            intercept + sum(coefficients[name] * state.outgoing
                for (name, state) in node.states))
    else
        @constraint(node.subproblem, bellman_function.variable <=
            intercept + sum(coefficients[name] * state.outgoing
                for (name, state) in node.states))
    end
    push!(bellman_function.cuts, Cut(intercept, coefficients))
    return
end
