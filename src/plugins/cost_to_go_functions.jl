#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

mutable struct Cut
    intercept::Float64
    coefficients::Dict{Symbol, Float64}
    non_dominated_count::Int
    constraint_ref::Union{Nothing, JuMP.ConstraintRef}
end

mutable struct SampledState
    state::Dict{Symbol, Float64}
    dominating_cut::Cut
    best_objective::Float64
end

struct LevelOneOracle
    cuts::Vector{Cut}
    states::Vector{SampledState}
    cuts_to_be_deleted::Vector{Cut}
    deletion_minimum::Int
    function LevelOneOracle(deletion_minimum)
        return new(Cut[], SampledState[], Cut[], deletion_minimum)
    end
end

struct ConvexApproximation
    θ::JuMP.VariableRef
    x′::Dict{Symbol, JuMP.VariableRef}
    cut_oracle::LevelOneOracle
    function ConvexApproximation(θ, x′; deletion_minimum)
        return new(θ, x′, LevelOneOracle(deletion_minimum))
    end
end

# Add the cut `V.θ ≥ θᵏ + ⟨πᵏ, x′ - xᵏ⟩`.
function _add_cut(V::ConvexApproximation, θᵏ, πᵏ, xᵏ)
    model = JuMP.owner_model(V.θ)
    for (key, x) in xᵏ
        θᵏ -= πᵏ[key] * xᵏ[key]
    end
    is_minimization = JuMP.objective_sense(model) == MOI.MIN_SENSE
    c_ref = if is_minimization
        @constraint(model, V.θ >= θᵏ + sum(πᵏ[i] * V.x′[i] for i in keys(V.x′)))
    else
        @constraint(model, V.θ <= θᵏ + sum(πᵏ[i] * V.x′[i] for i in keys(V.x′)))
    end
    cut = Cut(θᵏ, πᵏ, 1, c_ref)
    _level_one_update(V.cut_oracle, cut, xᵏ, is_minimization)
    return
end

# Internal function: calculate the height of `cut` evaluated at `state`.
function _eval_height(cut::Cut, state::Dict{Symbol, Float64})
    height = cut.intercept
    for (key, value) in cut.coefficients
        height += value * state[key]
    end
    return height
end

# Internal function: check if the candidate point dominates the incumbent.
function _dominates(candidate, incumbent, minimization::Bool)
    return minimization ? candidate > incumbent : candidate < incumbent
end

# Internal function: update the Level-One datastructures inside
# `cost_to_go`.
function _level_one_update(oracle::LevelOneOracle, cut::Cut,
                           state::Dict{Symbol, Float64}, is_minimization::Bool)
    sampled_state = SampledState(state, cut, _eval_height(cut, state))
    # Loop through previously sampled states and compare the height of the most
    # recent cut against the current best. If this cut is an improvement, store
    # this one instead.
    for old_state in oracle.states
        height = _eval_height(cut, old_state.state)
        if _dominates(height, old_state.best_objective, is_minimization)
            old_state.dominating_cut.non_dominated_count -= 1
            if old_state.dominating_cut.non_dominated_count <= 0
                push!(oracle.cuts_to_be_deleted, old_state.dominating_cut)
            end
            cut.non_dominated_count += 1
            old_state.dominating_cut = cut
            old_state.best_objective = height
        end
    end
    # Now loop through previously discovered cuts and compare their height at
    # the most recent sampled point in the state-space. If this cut is an
    # improvement, store this one instead. Note that we have to do this because
    # we might have previously thrown out a cut that is now relevant.
    for old_cut in oracle.cuts
        # If the constriant ref is not nothing, this cut is already in the
        # model, so it can't be better than the one we just found.
        if old_cut.constraint_ref !== nothing
            continue
        elseif !JuMP.is_valid(old_cut.constraint_ref)
            old_cut.constraint_ref = nothing
            continue
        end
        height = _eval_height(old_cut, state)
        if dominates(height, sampled_state.best_objective, is_minimization)
            sampled_state.dominating_cut.non_dominated_count -= 1
            if sampled_state.dominating_cut.non_dominated_count <= 0
                push!(oracle.cuts_to_be_deleted, sampled_state.dominating_cut)
            end
            old_cut.non_dominated_count += 1
            sampled_state.dominating_cut = old_cut
            sampled_state.best_objective = height
        end
    end
    push!(oracle.cuts, cut)
    push!(oracle.states, sampled_state)
    return
end

@enum(CutType, AVERAGE_CUT, MULTI_CUT)

struct CostToGoFunction <: AbstractBellmanFunction
    θ_global::ConvexApproximation
    θ_locals::Vector{ConvexApproximation}
    cut_type::CutType
end

cost_to_go_term(cost_to_go::CostToGoFunction) = cost_to_go.θ_global.θ

# Internal struct: this struct is just a cache for arguments until we can build
# an actual instance of the type T at a later point.
struct InstanceFactory{T}
    args
    kwargs
    InstanceFactory{T}(args...; kwargs...) where {T} = new{T}(args, kwargs)
end

function AverageCut(; lower_bound = -Inf, upper_bound = Inf,
                      deletion_minimum::Int = 1)
    return InstanceFactory{CostToGoFunction}(
        lower_bound = lower_bound, upper_bound = upper_bound,
        deletion_minimum = deletion_minimum)
end

function initialize_bellman_function(
        factory::InstanceFactory{CostToGoFunction}, model::PolicyGraph{T},
        node::Node{T}) where {T}
    lower_bound, upper_bound, deletion_minimum = -Inf, Inf, 0
    if length(factory.args) > 0
        error("Positional arguments $(factory.args) ignored in AverageCut.")
    end
    for (kw, value) in factory.kwargs
        if kw == :lower_bound
            lower_bound = value
        elseif kw == :upper_bound
            upper_bound = value
        elseif kw == :deletion_minimum
            deletion_minimum = value
        else
            error("Keyword $(kw) not recognised as argument to AverageCut.")
        end
    end
    if lower_bound == -Inf && upper_bound == Inf
        error("You must specify a finite bound on the cost-to-go term.")
    end
    if length(node.children) == 0
        lower_bound = upper_bound = 0.0
    end
    Θᴳ = @variable(node.subproblem)
    lower_bound > -Inf && JuMP.set_lower_bound(Θᴳ, lower_bound)
    upper_bound < Inf && JuMP.set_upper_bound(Θᴳ, upper_bound)
    x′ = Dict(key => var.out for (key, var) in node.states)
    return CostToGoFunction(
        ConvexApproximation(Θᴳ, x′; deletion_minimum=deletion_minimum),
        ConvexApproximation[],
        AVERAGE_CUT
    )
end

function refine_bellman_function(
            model::PolicyGraph{T},
            node::Node{T},
            cost_to_go::CostToGoFunction,
            risk_measure::AbstractRiskMeasure,
            outgoing_state::Dict{Symbol, Float64},
            dual_variables::Vector{Dict{Symbol, Float64}},
            noise_supports::Vector,
            nominal_probability::Vector{Float64},
            objective_realizations::Vector{Float64}) where {T}
    # Sanity checks.
    @assert length(dual_variables) == length(noise_supports) ==
        length(nominal_probability) == length(objective_realizations)
    # Preliminaries that are common to all cut types.
    risk_adjusted_probability = similar(nominal_probability)
    adjust_probability(
        risk_measure, risk_adjusted_probability, nominal_probability,
        noise_supports, objective_realizations,
        model.objective_sense == MOI.MIN_SENSE)
    # The meat of the function.
    if cost_to_go.cut_type == AVERAGE_CUT
        _add_average_cut(
            node,
            outgoing_state,
            risk_adjusted_probability,
            objective_realizations,
            dual_variables)
    else  # Add a multi-cut
        @assert cost_to_go.cut_type == MULTI_CUT
        _add_locals_if_necessary(cost_to_go, length(dual_variables))
        _add_multi_cut(
            node,
            outgoing_state,
            risk_adjusted_probability,
            objective_realizations,
            dual_variables)
    end
end

function _add_average_cut(node::Node, outgoing_state::Dict{Symbol, Float64},
                          risk_adjusted_probability::Vector{Float64},
                          objective_realizations::Vector{Float64},
                          duals::Vector{Dict{Symbol, Float64}})
    N = length(risk_adjusted_probability)
    @assert N == length(objective_realizations) == length(duals)
    # Calculate the expected intercept and dual variables with respect to the
    # risk-adjusted probability distributino.
    πᵏ = Dict(key => 0.0 for key in keys(outgoing_state))
    θᵏ = 0.0
    for i in 1:length(objective_realizations)
        p = risk_adjusted_probability[i]
        θᵏ += p * objective_realizations[i]
        for (key, dual) in duals[i]
            πᵏ[key] += p * dual
        end
    end
    # Now add the average-cut to the subproblem.
    _add_cut(node.bellman_function.θ_global, θᵏ, πᵏ, outgoing_state)
    return
end

function _add_multi_cut(node::Node, outgoing_state::Dict{Symbol, Float64},
                        risk_adjusted_probability::Vector{Float64},
                        objective_realizations::Vector{Float64},
                        duals::Vector{Dict{Symbol, Float64}})
    N = length(risk_adjusted_probability)
    @assert N == length(objective_realizations) == length(duals)
    cost_to_go = node.bellman_function
    for i in 1:length(duals)
        _add_cut(
            cost_to_go.θ_locals[i],
            objective_realizations[i], dual_variables[i], outgoing_state)
    end
    model = JuMP.owner_model(cost_to_go.θ_global)
    # TODO(odow): hash the risk_adjusted_probability and only add if it's a new
# probability distribution.
    if JuMP.objective_sense(model) == MOI.MIN_SENSE
        @constraint(model, cost_to_go.ϴ_global.θ >= sum(
            risk_adjusted_probability[i] * cost_to_go.θ_locals[i].θ
                for i in 1:length(risk_adjusted_probability)))
    else
        @constraint(model, cost_to_go.ϴ_global.θ <= sum(
            risk_adjusted_probability[i] * cost_to_go.θ_locals[i].θ
                for i in 1:length(risk_adjusted_probability)))
    end
    return
end

# If we are adding a multi-cut for the first time, then the local θ variables
# won't have been added.
# TODO(odow): a way to set different bounds for each variable in the multi-cut.
function _add_locals_if_necessary(cost_to_go::CostToGoFunction, N::Int)
    num_local_thetas = length(cost_to_go.θ_locals)
    if num_local_thetas == N
        # Do nothing. Already initialized.
    elseif num_local_thetas == 0
        Θᴳ = cost_to_go.θ_global.θ
        model = JuMP.owner_model(Θᴳ)
        for i in 1:N
            θ = @variable(model)
            if JuMP.has_lower_bound(Θᴳ)
                JuMP.set_lower_bound(θ, JuMP.lower_bound(Θᴳ))
            end
            if JuMP.has_upper_bound(cost_to_go.θ_global.θ)
                JuMP.set_upper_bound(θ, JuMP.upper_bound(Θᴳ))
            end
            push!(cost_to_go.θ_locals, ConvexApproximation(
                θ, Θᴳ.x′; deletion_minimum=Θᴳ.deletion_minimum))
        end
    else
        error("Expected $(N) local θ variables but there were $(num_local_thetas).")
    end
    return
end

"""
    write_cuts_to_file(model::PolicyGraph{T}, filename::String) where {T}

Write the cuts that form the policy in `model` to `filename` in JSON format.

See also [`SDDP.read_cuts_from_file`](@ref).
"""
function write_cuts_to_file(model::PolicyGraph{T}, filename::String) where {T}
    cuts = Dict{T, Vector{Dict{Symbol, Float64}}}()
    for (node_name, node) in model.nodes
        if node.objective_state !== nothing
            error("Unable to write cuts to file because it contains objective" *
                  " states.")
        end
        cuts[node_name] = Dict{String, Float64}[]
        for cut in node.bellman_function.θ_global.cut_oracle.cuts
            cut_dict = copy(cut.coefficients)
            cut_dict[:cut_intercept] = cut.intercept
            push!(cuts[node_name], cut_dict)
        end
    end
    open(filename, "w") do io
        write(io, JSON.json(cuts))
    end
    return
end

_node_name_parser(::Type{Int}, name::String) = parse(Int, name)
_node_name_parser(::Type{Symbol}, name::String) = Symbol(name)
function _node_name_parser(::Type{NTuple{N, Int}}, name::String) where {N}
    keys = parse.(Int, strip.(split(name[2:end-1], ",")))
    if length(keys) != N
        error("Unable to parse node called $(name). Expected $N elements.")
    end
    return tuple(keys...)
end

function _node_name_parser(T, name)
    error("Unable to read name $(name). Provide a custom parser to " *
          "`read_cuts_from_file` using the `node_name_parser` keyword.")
end

"""
    read_cuts_from_file(
        model::PolicyGraph{T}, filename::String;
        node_name_parser::Function = _node_name_parser) where {T}

Read cuts (saved using [`SDDP.write_cuts_to_file`](@ref)) from `filename` into
`model`.

Since `T` can be an arbitrary Julia type, the conversion to JSON is lossy. When
reading, `read_cuts_from_file` only supports `T=Int`, `T=NTuple{N, Int}`, and
`T=Symbol`. If you have manually created a policy graph with a different node
type `T`, provide a function `node_name_parser` with the signature
`node_name_parser(T, name::String)::T where {T}` that returns the name of each
node given the string name `name`.

See also [`SDDP.write_cuts_to_file`](@ref).
"""
function read_cuts_from_file(
        model::PolicyGraph{T}, filename::String;
        node_name_parser::Function = _node_name_parser) where {T}
    cuts = JSON.parsefile(filename, use_mmap=false)
    for (str_node_name, cut_list) in cuts
        node_name = node_name_parser(T, str_node_name)::T
        node = model[node_name]
        for json_cut in cut_list
            intercept = 0.0
            coefficients = Dict{Symbol, Float64}()
            for (name, coef) in json_cut
                if name == "cut_intercept"
                    intercept = coef
                else
                    coefficients[Symbol(name)] = coef
                end
            end
            _add_cut(
                node.bellman_function.θ_global,
                intercept,
                coefficients,
                Dict(key=>0.0 for key in keys(coefficients))
            )
        end
    end
    return
end
