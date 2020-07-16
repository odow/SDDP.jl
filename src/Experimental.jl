module Experimental

using SDDP

import JSON

function Base.write(io::IO, model::SDDP.PolicyGraph)
    edges = Dict{String, Any}[]
    _add_edges(edges, "$(model.root_node)", model.root_children)
    nodes = Dict{String, Any}()
    for (node_name, node) in model.nodes
        _add_edges(edges, "$(node_name)", node.children)
        _add_node_to_dict(nodes, node, node_name)
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
    )
    return Base.write(io, JSON.json(sof))
end

function write_to_file(model::SDDP.PolicyGraph, filename::String)
    return open(io -> Base.write(io, model), filename, "w")
end

function _add_edges(edges, from::String, children)
    for child in children
        push!(
            edges,
            Dict(
                "from" => from,
                "to" => "$(child.term)",
                "probability" => "$(child.probability)",
            )
        )
    end
end

function _add_node_to_dict(dest::Dict, node::SDDP.Node, node_name)
    realizations = Dict{Any, Any}[]
    random_variables = String[]
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
        else
            @assert length(setdiff(random_variables, keys(support))) == 0
        end
    end
    for x in random_variables
        unfix(variable_by_name(node.subproblem, x))
    end
    set_objective_function(node.subproblem, node.stage_objective)
    dest_model = MOI.FileFormats.Model(format = MOI.FileFormats.FORMAT_MOF)
    MOI.copy_to(dest_model, backend(node.subproblem))
    io = IOBuffer()
    Base.write(io, dest_model)
    seekstart(io)
    subproblem = JSON.parse(io)
    dest["$(node_name)"] = Dict(
        "state_variables" => Dict(
            "$(state_name)" => Dict(
                "in" => name(state.in), "out" => name(state.out)
            )
            for (state_name, state) in node.states
        ),
        "random_variables" => random_variables,
        "subproblem" => subproblem,
        "realizations" => realizations
    )
    return
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

function read_from_file(filename::String)
    data = JSON.parsefile(filename)
    graph = SDDP.Graph(data["root"]["name"])
    for node in keys(data["nodes"])
        SDDP.add_node(graph, node)
    end
    for edge in data["edges"]
        SDDP.add_edge(
            graph,
            edge["from"] => edge["to"],
            parse(Float64, edge["probability"])
        )
    end
    model = SDDP.PolicyGraph(
        graph,
        sense = :Min,
        lower_bound = -1e6,
    ) do sp, node_name
        _load_mof_model(sp, data, node_name)
        node = SDDP.get_node(sp)
        for (s, state) in data["nodes"][node_name]["state_variables"]
            node.states[Symbol(s)] = SDDP.State(
                variable_by_name(node.subproblem, state["in"]),
                variable_by_name(node.subproblem, state["out"]),
            )
        end
        Ω, P = Dict[], Float64[]
        for realization in data["nodes"][node_name]["realizations"]
            push!(P, realization["probability"])
            push!(Ω, realization["support"])
        end
        SDDP.parameterize(sp, Ω, P) do ω
            for (k, v) in ω
                fix(variable_by_name(sp, k), v)
            end
        end
        SDDP.set_stage_objective(sp, objective_function(sp))
    end
    for (k, v) in data["root"]["state_variables"]
        model.initial_root_state[Symbol(k)] = v["initial_value"]
    end
    return model
end

end
