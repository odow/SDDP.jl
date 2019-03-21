#  Copyright 2017-19, Oscar Dowson.
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
    function ConvexApproximation(θ, x′, deletion_minimum)
        return new(θ, x′, LevelOneOracle(deletion_minimum))
    end
end

# Add the cut `V.θ ≥ θᵏ + ⟨πᵏ, x′ - xᵏ⟩`.
function _add_cut(V::ConvexApproximation, θᵏ, πᵏ, xᵏ, μᵀy=JuMP.AffExpr(0.0); cut_selection::Bool=true)
    model = JuMP.owner_model(V.θ)
    for (key, x) in xᵏ
        θᵏ -= πᵏ[key] * xᵏ[key]
    end
    is_minimization = JuMP.objective_sense(model) == MOI.MIN_SENSE
    cut = Cut(θᵏ, πᵏ, 1, nothing)
    if μᵀy == JuMP.AffExpr(0.0) && cut_selection
        _level_one_update(V.cut_oracle, cut, xᵏ, is_minimization)
        _purge_cuts(V)
    end
    cut.constraint_ref = if is_minimization
        @constraint(model,
            V.θ + μᵀy >= θᵏ + sum(πᵏ[i] * V.x′[i] for i in keys(V.x′)))
    else
        @constraint(model,
            V.θ + μᵀy <= θᵏ + sum(πᵏ[i] * V.x′[i] for i in keys(V.x′)))
    end
    return
end

function _purge_cuts(V::ConvexApproximation)
    model = JuMP.owner_model(V.θ)
    if length(V.cut_oracle.cuts_to_be_deleted) >= V.cut_oracle.deletion_minimum
        for cut in V.cut_oracle.cuts_to_be_deleted
            JuMP.delete(model, cut.constraint_ref)
        end
        empty!(V.cut_oracle.cuts_to_be_deleted)
    end
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
# `bellman_function`.
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

# Internal struct: this struct is just a cache for arguments until we can build
# an actual instance of the type T at a later point.
struct InstanceFactory{T}
    args
    kwargs
    InstanceFactory{T}(args...; kwargs...) where {T} = new{T}(args, kwargs)
end

struct BellmanFunction <: AbstractBellmanFunction
    θ_global::ConvexApproximation
    θ_locals::Vector{ConvexApproximation}
    cut_type::CutType
end

function BellmanFunction(;
        lower_bound = -Inf, upper_bound = Inf, deletion_minimum::Int = 1,
        cut_type::CutType = AVERAGE_CUT)
    return InstanceFactory{BellmanFunction}(
        lower_bound = lower_bound, upper_bound = upper_bound,
        deletion_minimum = deletion_minimum, cut_type = cut_type)
end

bellman_term(bellman_function::BellmanFunction) = bellman_function.θ_global.θ

function initialize_bellman_function(
        factory::InstanceFactory{BellmanFunction}, model::PolicyGraph{T},
        node::Node{T}) where {T}
    lower_bound, upper_bound, deletion_minimum, cut_type = -Inf, Inf, 0, AVERAGE_CUT
    if length(factory.args) > 0
        error("Positional arguments $(factory.args) ignored in BellmanFunction.")
    end
    for (kw, value) in factory.kwargs
        if kw == :lower_bound
            lower_bound = value
        elseif kw == :upper_bound
            upper_bound = value
        elseif kw == :deletion_minimum
            deletion_minimum = value
        elseif kw == :cut_type
            cut_type = value
        else
            error("Keyword $(kw) not recognised as argument to BellmanFunction.")
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
    # Initialize bounds for the objective states. If objective_state==nothing,
    # this check will be skipped by dispatch.
    _add_initial_bounds(node.objective_state, Θᴳ)
    x′ = Dict(key => var.out for (key, var) in node.states)
    return BellmanFunction(ConvexApproximation(Θᴳ, x′, deletion_minimum),
                           ConvexApproximation[], cut_type)
end

# Internal function: helper used in _add_initial_bounds.
function _add_objective_state_constraint(
        θ::JuMP.VariableRef, y::NTuple{N, Float64},
        μ::NTuple{N, JuMP.VariableRef}) where {N}
    model = JuMP.owner_model(θ)
    lower_bound = JuMP.has_lower_bound(θ) ? JuMP.lower_bound(θ) : -Inf
    upper_bound = JuMP.has_upper_bound(θ) ? JuMP.upper_bound(θ) : Inf
    if lower_bound > -Inf
        @constraint(model, sum(y[i] * μ[i] for i in 1:N) + θ >= lower_bound)
    end
    if upper_bound < Inf
        @constraint(model, sum(y[i] * μ[i] for i in 1:N) + θ <= upper_bound)
    end
    if lower_bound ≈ upper_bound ≈ 0.0
        @constraint(model, [i=1:N], μ[i] == 0.0)
    end
    return
end

# Internal function: When created, θ has bounds of [-M, M], but, since we are
# adding these μ terms, we really want to bound <y, μ> + θ ∈ [-M, M]. We need to
# consider all possible values for `y`. Because the domain of `y` is
# rectangular, we want to add a constraint at each extreme point. This involves
# adding 2^N constraints where N = |μ|. This is only feasible for
# low-dimensional problems, e.g., N < 5.
_add_initial_bounds(obj_state::Nothing, θ) = nothing
function _add_initial_bounds(obj_state::ObjectiveState, θ)
    model = JuMP.owner_model(θ)
    if length(obj_state.μ) < 5
        for y in Base.product(zip(obj_state.lower_bound, obj_state.upper_bound)...)
            _add_objective_state_constraint(θ, y, obj_state.μ)
        end
    else
        _add_objective_state_constraint(θ, obj_state.lower_bound, obj_state.μ)
        _add_objective_state_constraint(θ, obj_state.upper_bound, obj_state.μ)
    end
end

function refine_bellman_function(
            model::PolicyGraph{T},
            node::Node{T},
            bellman_function::BellmanFunction,
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
    if bellman_function.cut_type == AVERAGE_CUT
        _add_average_cut(
            node,
            outgoing_state,
            risk_adjusted_probability,
            objective_realizations,
            dual_variables)
    else  # Add a multi-cut
        @assert bellman_function.cut_type == MULTI_CUT
        _add_locals_if_necessary(bellman_function, length(dual_variables))
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
    # Now add the average-cut to the subproblem. We include the objective-state
    # component μᵀy.
    _add_cut(node.bellman_function.θ_global, θᵏ, πᵏ, outgoing_state,
             get_objective_state_component(node))
    return
end

function _add_multi_cut(node::Node, outgoing_state::Dict{Symbol, Float64},
                        risk_adjusted_probability::Vector{Float64},
                        objective_realizations::Vector{Float64},
                        duals::Vector{Dict{Symbol, Float64}})
    N = length(risk_adjusted_probability)
    @assert N == length(objective_realizations) == length(duals)
    bellman_function = node.bellman_function
    for i in 1:length(duals)
        # Do not include the objective state component μᵀy.
        _add_cut(bellman_function.θ_locals[i], objective_realizations[i],
                 dual_variables[i], outgoing_state)
    end
    model = JuMP.owner_model(bellman_function.θ_global)
    μᵀy = get_objective_state_component(node)
    # TODO(odow): hash the risk_adjusted_probability and only add if it's a new
    # probability distribution.
    if JuMP.objective_sense(model) == MOI.MIN_SENSE
        @constraint(model, bellman_function.ϴ_global.θ + μᵀy >= sum(
            risk_adjusted_probability[i] * bellman_function.θ_locals[i].θ
                for i in 1:length(risk_adjusted_probability)))
    else
        @constraint(model, bellman_function.ϴ_global.θ + μᵀy <= sum(
            risk_adjusted_probability[i] * bellman_function.θ_locals[i].θ
                for i in 1:length(risk_adjusted_probability)))
    end
    return
end

# If we are adding a multi-cut for the first time, then the local θ variables
# won't have been added.
# TODO(odow): a way to set different bounds for each variable in the multi-cut.
function _add_locals_if_necessary(bellman_function::BellmanFunction, N::Int)
    num_local_thetas = length(bellman_function.θ_locals)
    if num_local_thetas == N
        # Do nothing. Already initialized.
    elseif num_local_thetas == 0
        Θᴳ = bellman_function.θ_global.θ
        model = JuMP.owner_model(Θᴳ)
        for i in 1:N
            θ = @variable(model)
            if JuMP.has_lower_bound(Θᴳ)
                JuMP.set_lower_bound(θ, JuMP.lower_bound(Θᴳ))
            end
            if JuMP.has_upper_bound(bellman_function.θ_global.θ)
                JuMP.set_upper_bound(θ, JuMP.upper_bound(Θᴳ))
            end
            push!(bellman_function.θ_locals, ConvexApproximation(
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
            error("Unable to write cuts to file because model contains " *
                  "objective states.")
        elseif length(node.bellman_function.θ_locals) > 0
            error("Unable to write cuts to file because model contains " *
                  "multi-cuts.")
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
            # So they cuts are written to file after they have been normalized
            # to `θᴳ ≥ [θᵏ - ⟨πᵏ, xᵏ⟩] + ⟨πᵏ, x′⟩`. Thus, we pass `xᵏ=0` so that
            # eveything works out okay.
            # Importantly, don't run cut selection when adding these cuts.
            _add_cut(
                node.bellman_function.θ_global,
                intercept,
                coefficients,
                Dict(key=>0.0 for key in keys(coefficients)),
                cut_selection = false
            )
        end
    end
    return
end
