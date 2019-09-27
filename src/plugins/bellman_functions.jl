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

mutable struct LevelOneOracle
    cuts::Vector{Cut}
    states::Vector{SampledState}
    cuts_to_be_deleted::Vector{Cut}
    deletion_minimum::Int
    function LevelOneOracle(deletion_minimum)
        return new(Cut[], SampledState[], Cut[], deletion_minimum)
    end
end

struct ConvexApproximation
    theta::JuMP.VariableRef
    states::Dict{Symbol, JuMP.VariableRef}
    cut_oracle::LevelOneOracle
    function ConvexApproximation(theta, states, deletion_minimum)
        return new(theta, states, LevelOneOracle(deletion_minimum))
    end
end

# Add the cut `V.θ ≥ θᵏ + ⟨πᵏ, x′ - xᵏ⟩`.
function _add_cut(
    V::ConvexApproximation,
    θᵏ::Float64,
    πᵏ::Dict{Symbol, Float64},
    xᵏ::Dict{Symbol, Float64},
    μᵀy = JuMP.AffExpr(0.0);
    cut_selection::Bool = true
)
    model = JuMP.owner_model(V.theta)
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
        @constraint(
            model,
            V.theta + μᵀy >= θᵏ + sum(πᵏ[i] * x for (i, x) in V.states)
        )
    else
        @constraint(
            model,
            V.theta + μᵀy <= θᵏ + sum(πᵏ[i] * x for (i, x) in V.states)
        )
    end
    return
end

function _purge_cuts(V::ConvexApproximation)
    model = JuMP.owner_model(V.theta)
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
function _level_one_update(
    oracle::LevelOneOracle,
    cut::Cut,
    state::Dict{Symbol, Float64},
    is_minimization::Bool
)
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

@enum(CutType, SINGLE_CUT, MULTI_CUT)

# Internal struct: this struct is just a cache for arguments until we can build
# an actual instance of the type T at a later point.
struct InstanceFactory{T}
    args
    kwargs
    InstanceFactory{T}(args...; kwargs...) where {T} = new{T}(args, kwargs)
end

mutable struct BellmanFunction <: AbstractBellmanFunction
    global_theta::ConvexApproximation
    local_thetas::Vector{ConvexApproximation}
    cut_type::CutType
    # Cuts defining the dual representation of the risk measure.
    risk_set_cuts::Set{Vector{Float64}}
end

"""
    BellmanFunction(;
        lower_bound = -Inf, upper_bound = Inf, deletion_minimum::Int = 1,
        cut_type::CutType = MULTI_CUT)
"""
function BellmanFunction(
    ;
    lower_bound = -Inf,
    upper_bound = Inf,
    deletion_minimum::Int = 1,
    cut_type::CutType = MULTI_CUT
)
    return InstanceFactory{BellmanFunction}(
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        deletion_minimum = deletion_minimum,
        cut_type = cut_type
    )
end

function bellman_term(bellman_function::BellmanFunction)
    return bellman_function.global_theta.theta
end

function initialize_bellman_function(
    factory::InstanceFactory{BellmanFunction},
    model::PolicyGraph{T},
    node::Node{T}
) where {T}
    lower_bound, upper_bound, deletion_minimum, cut_type = -Inf, Inf, 0, SINGLE_CUT
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
    return BellmanFunction(
        ConvexApproximation(Θᴳ, x′, deletion_minimum),
        ConvexApproximation[],
        cut_type,
        Set{Vector{Float64}}()
    )
end

# Internal function: helper used in _add_initial_bounds.
function _add_objective_state_constraint(
    theta::JuMP.VariableRef,
    y::NTuple{N, Float64},
    μ::NTuple{N, JuMP.VariableRef}
) where {N}
    model = JuMP.owner_model(theta)
    lower_bound = JuMP.has_lower_bound(theta) ? JuMP.lower_bound(theta) : -Inf
    upper_bound = JuMP.has_upper_bound(theta) ? JuMP.upper_bound(theta) : Inf
    if lower_bound > -Inf
        @constraint(model, sum(y[i] * μ[i] for i in 1:N) + theta >= lower_bound)
    end
    if upper_bound < Inf
        @constraint(model, sum(y[i] * μ[i] for i in 1:N) + theta <= upper_bound)
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
_add_initial_bounds(obj_state::Nothing, theta) = nothing
function _add_initial_bounds(obj_state::ObjectiveState, theta)
    model = JuMP.owner_model(theta)
    if length(obj_state.μ) < 5
        for y in Base.product(zip(obj_state.lower_bound, obj_state.upper_bound)...)
            _add_objective_state_constraint(theta, y, obj_state.μ)
        end
    else
        _add_objective_state_constraint(theta, obj_state.lower_bound, obj_state.μ)
        _add_objective_state_constraint(theta, obj_state.upper_bound, obj_state.μ)
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
    objective_realizations::Vector{Float64}
) where {T}
    # Sanity checks.
    @assert length(dual_variables) == length(noise_supports) ==
        length(nominal_probability) == length(objective_realizations)
    # Preliminaries that are common to all cut types.
    risk_adjusted_probability = similar(nominal_probability)
    adjust_probability(
        risk_measure,
        risk_adjusted_probability,
        nominal_probability,
        noise_supports,
        objective_realizations,
        model.objective_sense == MOI.MIN_SENSE
    )
    # The meat of the function.
    if bellman_function.cut_type == SINGLE_CUT
        _add_average_cut(
            node,
            outgoing_state,
            risk_adjusted_probability,
            objective_realizations,
            dual_variables
        )
    else  # Add a multi-cut
        @assert bellman_function.cut_type == MULTI_CUT
        _add_locals_if_necessary(bellman_function, length(dual_variables))
        _add_multi_cut(
            node,
            outgoing_state,
            risk_adjusted_probability,
            objective_realizations,
            dual_variables
        )
    end
end

function _add_average_cut(
    node::Node,
    outgoing_state::Dict{Symbol, Float64},
    risk_adjusted_probability::Vector{Float64},
    objective_realizations::Vector{Float64},
    dual_variables::Vector{Dict{Symbol, Float64}}
)
    N = length(risk_adjusted_probability)
    @assert N == length(objective_realizations) == length(dual_variables)
    # Calculate the expected intercept and dual variables with respect to the
    # risk-adjusted probability distributino.
    πᵏ = Dict(key => 0.0 for key in keys(outgoing_state))
    θᵏ = 0.0
    for i in 1:length(objective_realizations)
        p = risk_adjusted_probability[i]
        θᵏ += p * objective_realizations[i]
        for (key, dual) in dual_variables[i]
            πᵏ[key] += p * dual
        end
    end
    # Now add the average-cut to the subproblem. We include the objective-state
    # component μᵀy and the belief state (if it exists).
    μᵀy = get_objective_state_component(node)
    JuMP.add_to_expression!(μᵀy, get_belief_state_component(node))
    _add_cut(node.bellman_function.global_theta, θᵏ, πᵏ, outgoing_state, μᵀy)
    return
end

function _add_multi_cut(
    node::Node,
    outgoing_state::Dict{Symbol, Float64},
    risk_adjusted_probability::Vector{Float64},
    objective_realizations::Vector{Float64},
    dual_variables::Vector{Dict{Symbol, Float64}}
)
    N = length(risk_adjusted_probability)
    @assert N == length(objective_realizations) == length(dual_variables)
    bellman_function = node.bellman_function
    μᵀy = get_objective_state_component(node)
    JuMP.add_to_expression!(μᵀy, get_belief_state_component(node))
    for i in 1:length(dual_variables)
        _add_cut(
            bellman_function.local_thetas[i],
            objective_realizations[i],
            dual_variables[i],
            outgoing_state,
            μᵀy
        )
    end
    model = JuMP.owner_model(bellman_function.global_theta.theta)
    cut_expr = @expression(
        model,
        sum(
            risk_adjusted_probability[i] * bellman_function.local_thetas[i].theta
            for i in 1:N
        ) - (1 - sum(risk_adjusted_probability)) * μᵀy
    )
    # TODO(odow): should we use `cut_expr` instead?
    ξ = copy(risk_adjusted_probability)
    if !(ξ in bellman_function.risk_set_cuts) || μᵀy != JuMP.AffExpr(0.0)
        push!(bellman_function.risk_set_cuts, ξ)
        if JuMP.objective_sense(model) == MOI.MIN_SENSE
            @constraint(model, bellman_function.global_theta.theta >= cut_expr)
        else
            @constraint(model, bellman_function.global_theta.theta <= cut_expr)
        end
    end
    return
end

# If we are adding a multi-cut for the first time, then the local θ variables
# won't have been added.
# TODO(odow): a way to set different bounds for each variable in the multi-cut.
function _add_locals_if_necessary(bellman_function::BellmanFunction, N::Int)
    num_local_thetas = length(bellman_function.local_thetas)
    if num_local_thetas == N
        # Do nothing. Already initialized.
    elseif num_local_thetas == 0
        global_theta = bellman_function.global_theta
        model = JuMP.owner_model(global_theta.theta)
        local_thetas = @variable(model, [1:N])
        if JuMP.has_lower_bound(global_theta.theta)
            JuMP.set_lower_bound.(local_thetas, JuMP.lower_bound(global_theta.theta))
        end
        if JuMP.has_upper_bound(global_theta.theta)
            JuMP.set_upper_bound.(local_thetas, JuMP.upper_bound(global_theta.theta))
        end
        for local_theta in local_thetas
            push!(
                bellman_function.local_thetas,
                ConvexApproximation(
                    local_theta, global_theta.states,
                    global_theta.cut_oracle.deletion_minimum
                )
            )
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
    cuts = Dict{String, Any}[]
    for (node_name, node) in model.nodes
        if node.objective_state !== nothing
            error("Unable to write cuts to file because model contains " *
                  "objective states.")
        end
        node_cuts = Dict(
            "node" => string(node_name),
            "single_cuts" => Dict{String, Any}[],
            "multi_cuts" => Dict{String, Any}[],
            "risk_set_cuts" => Vector{Float64}[]
        )
        for cut in node.bellman_function.global_theta.cut_oracle.cuts
            push!(node_cuts["single_cuts"], Dict(
                "intercept" => cut.intercept,
                "coefficients" => copy(cut.coefficients)
            ))
        end
        for (i, theta) in enumerate(node.bellman_function.local_thetas)
            for cut in theta.cut_oracle.cuts
                push!(node_cuts["multi_cuts"], Dict(
                    "realization" => i,
                    "intercept" => cut.intercept,
                    "coefficients" => copy(cut.coefficients)
                ))
            end
        end
        for p in node.bellman_function.risk_set_cuts
            push!(node_cuts["risk_set_cuts"], p)
        end
        push!(cuts, node_cuts)
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
    model::PolicyGraph{T},
    filename::String;
    node_name_parser::Function = _node_name_parser
) where {T}
    # So the cuts are written to file after they have been normalized
    # to `θᴳ ≥ [θᵏ - ⟨πᵏ, xᵏ⟩] + ⟨πᵏ, x′⟩`. Thus, we pass `xᵏ=0` so that
    # eveything works out okay.
    # Importantly, don't run cut selection when adding these cuts.
    cuts = JSON.parsefile(filename, use_mmap=false)
    for node_cuts in cuts
        node_name = node_name_parser(T, node_cuts["node"])::T
        node = model[node_name]
        bf = node.bellman_function
        # Loop through and add the single-cuts.
        for json_cut in node_cuts["single_cuts"]
            coefficients = Dict{Symbol, Float64}(
                Symbol(k) => v for (k, v) in json_cut["coefficients"])
            _add_cut(
                bf.global_theta,
                json_cut["intercept"],
                coefficients,
                Dict(key=>0.0 for key in keys(coefficients)),
                cut_selection = false
            )
        end
        # Loop through and add the multi-cuts. There are two parts:
        #  (i) the cuts w.r.t. the state variable x
        # (ii) the cuts that define the risk set
        # There is one additional complication: if these cuts are being read
        # into a new model, the local theta variables may not exist yet.
        if length(node_cuts["risk_set_cuts"]) > 0
            _add_locals_if_necessary(
                bf, length(first(node_cuts["risk_set_cuts"])))
        end
        for json_cut in node_cuts["multi_cuts"]
            coefficients = Dict{Symbol, Float64}(
                Symbol(k) => v for (k, v) in json_cut["coefficients"])
            _add_cut(
                bf.local_thetas[json_cut["realization"]],
                json_cut["intercept"],
                coefficients,
                Dict(key=>0.0 for key in keys(coefficients)),
                cut_selection = false
            )
        end
        # Here is part (ii): adding the constraints that define the risk-set
        # representation of the risk measure.
        for json_cut in node_cuts["risk_set_cuts"]
            expr = @expression(node.subproblem,
                bf.global_theta.theta - sum(
                    p * V.theta for (p, V) in zip(json_cut, bf.local_thetas)
                )
            )
            if JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE
                @constraint(node.subproblem, expr >= 0)
            else
                @constraint(node.subproblem, expr <= 0)
            end
        end
    end
    return
end
