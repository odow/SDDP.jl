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
    random_variables, realizations = String, Dict{String, Any}[]
    for noise in node.noise_terms
        support = Dict{String, Float64}()
        node.parameterize(noise.term)
        for x in all_variables(node.subproblem)
            if is_fixed(x)
                support[name(x)] = fix_value(x)
            end
        end
        push!(
            realizations,
            Dict("probability" => noise.probability, "support" => support)
        )
        if length(realizations) == 1
            random_variables = collect(keys(support))
        end
        @assert length(setdiff(random_variables, keys(support))) == 0
    end
    for x in random_variables
        unfix(variable_by_name(node.subproblem, x))
    end
    added_variables = _reformulate_objective_uncertainty(
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
    if length(added_variables) > 0
        # Undo the reformulation in user-space.
        delete(node.subproblem, added_variables)
        set_objective_function(node.subproblem, node.stage_objective)
    end
    return
end

function _reformulate_objective_uncertainty(
    node::Node, realizations, random_variables
)
    # TODO(odow): handle quadratic objectives.
    # First loop:
    #   - build different objectives
    #   - detect ones that change
    objectives = AffExpr[]
    constant_changes, coefficient_changes = false, Set{VariableRef}()
    for noise in node.noise_terms
        node.parameterize(noise.term)
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
    objective = copy(node.stage_objective)
    # Reformulate a changing objective constant.
    if constant_changes
        constant_term = @variable(node.subproblem)
        push!(added_variables, constant_term)
        set_name(constant_term, "_SDDPjl_random_objective_constant_")
        push!(random_variables, "_SDDPjl_random_objective_constant_")
        for (r, o) in zip(realizations, objectives)
            r["support"]["_SDDPjl_random_objective_constant_"] = o.constant
        end
        objective.constant = 0.0
        objective.terms[constant_term] = 1.0
    end
    # Reformulate changing objective coefficients.
    if length(coefficient_changes) > 0
        objective = convert(QuadExpr, objective)
        for x in coefficient_changes
            coef_term = @variable(node.subproblem)
            push!(added_variables, coef_term)
            coef_name = "_SDDPjl_random_objective_$(name(x))_"
            set_name(coef_term, coef_name)
            push!(random_variables, coef_name)
            for (r, o) in zip(realizations, objectives)
                r["support"][coef_name] = get(o.terms, x, 0.0)
            end
            delete!.(Ref(objective.aff.terms), x)
            add_to_expression!(objective, 1.0, coef_term, x)
        end
    end
    # Set the objective function to be written out.
    set_objective_function(node.subproblem, objective)
    return added_variables
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
            sp, String.(data["nodes"][node_name]["random_variables"])
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
