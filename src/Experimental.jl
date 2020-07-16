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
    set_objective_function(node.subproblem, node.stage_objective)
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
    return
end

function _subproblem_to_dict(subproblem::Model)
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
        parameterize(sp, Ω, P) do ω
            for (k, v) in ω
                fix(variable_by_name(sp, k), v)
            end
        end
        if objective_sense(sp) == model_sense
            set_stage_objective(sp, objective_function(sp))
        else
            @warn(
                "Flipping the objective sense of node $(node_name) so that " *
                "it matches the majority of the subproblems."
            )
            set_stage_objective(sp, -objective_function(sp))
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
