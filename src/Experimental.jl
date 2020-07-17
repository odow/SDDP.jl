#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    Base.write(io::IO, model::PolicyGraph)

WARNING: THIS FUNCTION IS EXPERIMENTAL. THINGS MAY CHANGE BETWEEN COMMITS. YOU
SHOULD NOT RELY ON THIS FUNCTIONALITY AS A LONG-TERM FILE FORMAT (YET).

Write `model` to `io` in the StochOptFormat file format.
"""
function Base.write(io::IO, model::PolicyGraph)
    edges = Dict{String, Any}[]
    _add_edges(edges, "$(model.root_node)", model.root_children)
    nodes = Dict{String, Any}()
    for (node_name, node) in model.nodes
        _add_edges(edges, "$(node_name)", node.children)
        _add_node_to_dict(nodes, node, "$(node_name)")
    end
    sof = Dict(
        "description" => "A problem exported from SDDP.jl",
        "version" => Dict("major" => 0, "minor" => 1),
        "root" => Dict(
            "name" => "$(model.root_node)",
            "state_variables" => Dict(
                "$(k)" => Dict("initial_value" => v)
                for (k, v) in model.initial_root_state
            )
        ),
        "nodes" => nodes,
        "edges" => edges,
        # TODO(odow): generate `test_scenarios`.
        "test_scenarios" => Any[]
    )
    return Base.write(io, JSON.json(sof))
end

function _add_edges(
    edges::Vector{Dict{String, Any}}, from::String, children::Vector{<:Noise}
)
    for child in children
        push!(
            edges,
            Dict(
                "from" => from,
                "to" => "$(child.term)",
                "probability" => child.probability,
            )
        )
    end
end

function _add_node_to_dict(dest::Dict, node::Node, node_name::String)
    random_variables = String[]
    realizations = Dict{String, Any}[
        Dict{String, Any}(
            "probability" => noise.probability,
            "support" => Dict{String, Float64}()
        ) for noise in node.noise_terms
    ]
    undo_reformulation = _reformulate_uncertainty(
        node, realizations, random_variables
    )
    dest[node_name] = Dict(
        "state_variables" => Dict(
            "$(state_name)" => Dict(
                "in" => name(state.in), "out" => name(state.out)
            )
            for (state_name, state) in node.states
        ),
        "random_variables" => random_variables,
        "subproblem" => _subproblem_to_dict(node.subproblem),
        "realizations" => realizations,
    )
    undo_reformulation()
    return
end

"""
    _reformulate_uncertainty(
        node::Node, realizations, random_variables
    )

Convert any lower and upper bounds than depend on the uncertainty into linear
constraints with a random variable.

Fixed variables are recorded as random variables, but no transformation is done.

Given an affine stageobjective with stagewise independent uncertainty,
reformulate into a quadratic stage objective by replacing the random
coefficients with random decision variables.

Return a function that undoes the reformulation when called with no arguments.
"""
function _reformulate_uncertainty(
    node::Node, realizations::Vector, random_variables::Vector{String}
)
    # Storage for things that are changing.
    bounds = Dict{VariableRef, Any}[]
    lb = Set{VariableRef}()
    ub = Set{VariableRef}()
    fx = Set{VariableRef}()
    objectives = AffExpr[]
    constant_changes = false
    coefficient_changes = Set{VariableRef}()

    # Loop 1: Collect terms that are changing
    for noise in node.noise_terms
        node.parameterize(noise.term)

        # Collect changing variable bounds.
        bound = Dict{VariableRef, Any}()
        for x in all_variables(node.subproblem)
            l, u, f = -Inf, Inf, 0.0
            if has_lower_bound(x)
                l = lower_bound(x)
            end
            if has_upper_bound(x)
                u = upper_bound(x)
            end
            if is_fixed(x)
                f = fix_value(x)
            end
            if length(bounds) >= 1
                if bounds[1][x].l != l
                    push!(lb, x)
                end
                if bounds[1][x].u != u
                    push!(ub, x)
                end
                if bounds[1][x].f != f
                    push!(fx, x)
                end
            end
            bound[x] = (l = l, u = u, f = f)
        end
        push!(bounds, bound)

        # Collect changing objective terms.
        push!(objectives, convert(AffExpr, node.stage_objective))
        if length(objectives) > 1
            obj = objectives[end]
            if obj.constant != objectives[1].constant
                constant_changes = true
            end
            for k in _dict_diff_keys(objectives[1].terms, obj.terms)
                push!(coefficient_changes, k)
            end
        end
    end

    added_variables = VariableRef[]
    added_constraints = ConstraintRef[]
    objective = copy(node.stage_objective)

    # Reformulate a changing objective constant.
    if constant_changes
        new_name = "_SDDPjl_random_objective_constant_"
        y = _add_new_random_variable(
            node, new_name, random_variables, added_variables
        )
        for (r, o) in zip(realizations, objectives)
            r["support"][new_name] = o.constant
        end
        objective.constant = 0.0
        objective.terms[y] = 1.0
    end

    # Reformulate changing objective coefficients.
    if length(coefficient_changes) > 0
        objective = convert(QuadExpr, objective)
        for x in coefficient_changes
            new_name = "_SDDPjl_random_objective_$(name(x))_"
            y = _add_new_random_variable(
                node, new_name, random_variables, added_variables
            )
            for (r, o) in zip(realizations, objectives)
                r["support"][new_name] = get(o.terms, x, 0.0)
            end
            delete!.(Ref(objective.aff.terms), x)
            add_to_expression!(objective, 1.0, y, x)
        end
    end

    # Reformulate fixed variables.
    for x in fx
        for (realization, bound) in zip(realizations, bounds)
            realization["support"][name(x)] = bound[x].f
        end
        push!(random_variables, name(x))
        unfix(x)
    end

    # Reformulate lower bounded variables.
    for x in lb
        new_name = "_SDDPjl_lower_bound_$(name(x))_"
        y = _add_new_random_variable(
            node, new_name, random_variables, added_variables
        )
        c = @constraint(node.subproblem, x >= y)
        push!(added_constraints, c)
        delete_lower_bound(x)
        for (realization, bound) in zip(realizations, bounds)
            realization["support"][new_name] = bound[x].l
        end
    end

    # Reformulate upper bounded variables.
    for x in ub
        new_name = "_SDDPjl_upper_bound_$(name(x))_"
        y = _add_new_random_variable(
            node, new_name, random_variables, added_variables
        )
        c = @constraint(node.subproblem, x <= y)
        push!(added_constraints, c)
        delete_upper_bound(x)
        for (realization, bound) in zip(realizations, bounds)
            realization["support"][new_name] = bound[x].u
        end
    end

    # Set the objective function to be written out.
    set_objective_function(node.subproblem, objective)

    return () -> begin
        delete(node.subproblem, added_variables)
        delete.(Ref(node.subproblem), added_constraints)
        set_objective_function(node.subproblem, node.stage_objective)
        return
    end
end

function _add_new_random_variable(
    node, new_name, random_variables, added_variables
)
    y = @variable(node.subproblem, base_name = new_name)
    push!(added_variables, y)
    push!(random_variables, new_name)
    return y
end


function _dict_diff_keys(
    x::AbstractDict{K, V}, y::AbstractDict{K, V}
) where {K, V}
    diff = Set{K}()
    for (k, v) in x
        if haskey(y, k)
            if v != y[k]
                push!(diff, k)
            end
        else
            push!(diff, k)
        end
    end
    for k in keys(y)
        if !haskey(x, k)
            push!(diff, k)
        end
    end
    return diff
end

function _subproblem_to_dict(subproblem::JuMP.Model)
    dest_model = MOI.FileFormats.Model(format = MOI.FileFormats.FORMAT_MOF)
    MOI.copy_to(dest_model, backend(subproblem))
    io = IOBuffer()
    Base.write(io, dest_model)
    seekstart(io)
    return JSON.parse(io; dicttype = Dict{String, Any})
end

function _load_mof_model(sp::JuMP.Model, data::Dict, node::String)
    model = MOI.FileFormats.Model(format = MOI.FileFormats.FORMAT_MOF)
    io = IOBuffer()
    Base.write(io, JSON.json(data["nodes"][node]["subproblem"]))
    seekstart(io)
    MOI.read!(io, model)
    MOI.copy_to(sp, model)
    return
end

"""
    Base.read(io::IO, ::Type{PolicyGraph}; bound::Float64 = 1e6)

WARNING: THIS FUNCTION IS EXPERIMENTAL. THINGS MAY CHANGE BETWEEN COMMITS. YOU
SHOULD NOT RELY ON THIS FUNCTIONALITY AS A LONG-TERM FILE FORMAT (YET).

Return a `PolicyGraph` object read from `io` in the StochOptFormat file
format.
"""
function Base.read(io::IO, ::Type{PolicyGraph}; bound::Float64 = 1e6)
    data = JSON.parse(io; dicttype = Dict{String, Any})
    graph = Graph(data["root"]["name"])
    for (node_name, _) in data["nodes"]
        add_node(graph, node_name)
    end
    for edge in data["edges"]
        add_edge(graph, edge["from"] => edge["to"], edge["probability"])
    end
    proportion_min = sum(
        node["subproblem"]["objective"]["sense"] == "min"
        for (_, node) in data["nodes"]
    ) / length(data["nodes"])
    model_sense = proportion_min >= 0.5 ? MOI.MIN_SENSE : MOI.MAX_SENSE
    function subproblem_builder(sp::Model, node_name::String)
        _load_mof_model(sp, data, "$(node_name)")
        node = get_node(sp)
        for (s, state) in data["nodes"][node_name]["state_variables"]
            node.states[Symbol(s)] = State(
                variable_by_name(node.subproblem, state["in"]),
                variable_by_name(node.subproblem, state["out"]),
            )
        end
        Ω, P = Dict[], Float64[]
        for realization in data["nodes"][node_name]["realizations"]
            push!(P, realization["probability"])
            push!(Ω, realization["support"])
        end
        if objective_sense(sp) != model_sense
            @warn(
                "Flipping the objective sense of node $(node_name) so that " *
                "it matches the majority of the subproblems."
            )
        end
        obj_sgn = objective_sense(sp) == model_sense ? 1 : -1
        objective_coefficients, objf = _convert_objective_function(
            sp,
            convert(
                Vector{String},
                data["nodes"][node_name]["random_variables"]
            )
        )
        parameterize(sp, Ω, P) do ω
            for (k, v) in ω
                x = get(objective_coefficients, k, nothing)
                if x !== nothing
                    if objf isa AffExpr
                        objf.terms[x.var] = x.aff + v * x.coef
                    else
                        objf.aff.terms[x.var] = x.aff + v * x.coef
                    end
                end
                fix(variable_by_name(sp, k), v)
            end
            @stageobjective(sp, obj_sgn * objf)
        end
    end
    model = if model_sense == MOI.MIN_SENSE
        PolicyGraph(
            subproblem_builder, graph; sense = :Min, lower_bound = -abs(bound)
        )
    else
        PolicyGraph(
            subproblem_builder, graph; sense = :Max, upper_bound = abs(bound)
        )
    end
    for (k, v) in data["root"]["state_variables"]
        model.initial_root_state[Symbol(k)] = v["initial_value"]
    end
    return model
end

function _convert_objective_function(sp::Model, rvs::Vector{String})
    return _convert_objective_function(sp, rvs, objective_function(sp))
end

function _convert_objective_function(sp::Model, ::Vector{String}, objf)
    return Dict{String, Any}(), objf
end

function _convert_objective_function(
    sp::Model, rvs::Vector{String}, objf::QuadExpr
)
    terms = Dict{String, Any}()
    aff_obj = copy(objf.aff)
    quad_terms = empty(copy(objf.terms))
    for (k, v) in objf.terms
        a, b = name(k.a), name(k.b)
        if a in rvs && b in rvs
            error(
                "Please open an issue to support random * random in objective."
            )
        elseif a in rvs
            terms[a] = (var = k.b, coef = v, aff = get(aff_obj.terms, a, 0.0))
        elseif b in rvs
            terms[a] = (var = k.a, coef = v, aff = get(aff_obj.terms, b, 0.0))
        else
            quad_terms[k] = v
        end
    end
    if length(terms) == length(objf.terms)
        return terms, aff_obj
    else
        return terms, QuadExpr(aff_obj, quad_terms)
    end
end

"""
    write_to_file(model::PolicyGraph, filename::String)

WARNING: THIS FUNCTION IS EXPERIMENTAL. THINGS MAY CHANGE BETWEEN COMMITS. YOU
SHOULD NOT RELY ON THIS FUNCTIONALITY AS A LONG-TERM FILE FORMAT (YET).

Write `model` to `filename` in the StochOptFormat file format.
"""
function write_to_file(model::PolicyGraph, filename::String)
    return open(io -> Base.write(io, model), filename, "w")
end

"""
    read_from_file(filename::String; kwargs...)::PolicyGraph{String}

WARNING: THIS FUNCTION IS EXPERIMENTAL. THINGS MAY CHANGE BETWEEN COMMITS. YOU
SHOULD NOT RELY ON THIS FUNCTIONALITY AS A LONG-TERM FILE FORMAT (YET).

Return a `PolicyGraph` object read from `filename` in the StochOptFormat file
format.

See [`Base.read(::IO, ::Type{PolicyGraph})`](@ref) for information on the
keyword arguments that can be provided.
"""
function read_from_file(filename::String; kwargs...)
    return open(io -> Base.read(io, PolicyGraph; kwargs...), filename, "r")
end
