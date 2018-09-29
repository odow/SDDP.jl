abstract type AbstractBellmanFunction end

"""
    refine_bellman_function(graph::PolicyGraph{T},
                            node::Node{T},
                            bellman_function::AbstractBellmanFunction,
                            state::Dict{Symbol, Float64},
                            dual_variables::Vector{Dict{Symbol, Float64}},
                            noise_supports::Vector{<:Noise},
                            original_probability::Vector{Float64},
                            objective_realizations::Vector{Float64}
                                ) where T
"""
function refine_bellman_function end

function initialize(::Type{<:AbstractBellmanFunction}, graph::PolicyGraph{T},
                    node::Node{T}) where T
    error("initialize")
end

"""
    bellman_term(::AbstractBellmanFunction)

Return a JuMP expression representing the Bellman function.
"""
function bellman_term end

# ============================== SDDP.AverageCut ===============================

struct AverageCut <: AbstractBellmanFunction
    variable::JuMP.Variable
end

function initialize(
        ::Type{AverageCut}, graph::PolicyGraph{T}, node::Node{T}) where T
    node.bellman_function = AverageCut(@variable(node.subproblem))
    return
end

bellman_term(bellman::AverageCut) = bellman.variable

function refine_bellman_function(graph::PolicyGraph{T},
                                 node::Node{T},
                                 bellman_function::AbstractBellmanFunction,
                                 outgoing_state::Dict{Symbol, Float64},
                                 dual_variables::Vector{Dict{Symbol, Float64}},
                                 noise_supports::Vector{<:Noise},
                                 original_probability::Vector{Float64},
                                 objective_realizations::Vector{Float64}
                                     ) where T
    risk_adjusted_probability = similar(original_probability)
    adjust_probability(Expectation(),
                       risk_adjusted_probability,
                       original_probability,
                       noise_supports,
                       objective_realizations,
                       JuMP.objective_sense(node.subproblem) == :Min)
    @constraint(node.subproblem,
        sum(prob * (objective + sum(
                dual[x] * (state.outgoing - outgoing_state[name])
                    for (name, state) in node.states)
            for (objective, dual, prob) in zip(
                objective_realizations, dual_variables, risk_adjusted_probability)
        )
    )
end
