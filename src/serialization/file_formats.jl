#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP
using JSON
const PolicyGraph = SDDP.PolicyGraph

function _sof_build_root(model::PolicyGraph)
    root = Dict{String, Any}(
        "name" => "$(model.root_node)",
        "states" => Dict{String, Any}("$(k)" => v for (k, v) in model.initial_root_state)
    )
    return root
end

function _sof_build_nodes(model::PolicyGraph)
    nodes = Dict{String, Any}()
    for (k, node) in model.nodes
        io = IOBuffer()
        write(io, node.subproblem; format = MOI.FileFormats.FORMAT_MOF)
        seekstart(io)
        nodes["$(k)"] = Dict{String, Any}(
            "subproblem" => JSON.parse(io),
            "states" => Dict{String, Any}(
                "$(name)" => Dict{String, Any}(
                    "in" => "$(state.in)",
                    "out" => "$(state.out)"
                )
                for (name, state) in node.states
            ),
            "noise_terms" => Any[]
        )
    end
    return nodes
end

function _sof_build_edges(model::PolicyGraph)
    edges = Dict{String, Any}[
        Dict{String, Any}(
            "from" => "$(model.root_node)",
            "to" => "$(child.term)",
            "probability" => child.probability
        )
        for child in model.root_children
    ]
    for (k, node) in model.nodes
        for child in node.children
            push!(
                edges,
                Dict{String,Any}(
                    "from" => "$(k)",
                    "to" => "$(child.term)",
                    "probability" => child.probability
                )
            )
        end
    end
    return edges
end

function _sof_belief_partition(out::Dict{String, Any}, model::PolicyGraph)
    if model.belief_partition !== nothing
        out["belief_partition"] = model.belief_partition
    end
    return
end

function write_to_file(model::PolicyGraph, filename::String)
    out = Dict{String, Any}(
        "root" => _sof_build_root(model),
        "nodes" => _sof_build_nodes(model),
        "edges" => _sof_build_edges(model)
    )
    _sof_belief_partition(out, model)
    open(filename, "w") do io
        write(io, JSON.json(out, 2))
    end
    return
end

graph = SDDP.Graph(
    :root_node,
    [:week],
    [(:root_node => :week, 1.0), (:week => :week, 0.9)]
)
SDDP.add_ambiguity_set(graph, [:week])
model = SDDP.PolicyGraph(
    graph,
    lower_bound = 0
) do subproblem, node
    @variable(subproblem, x[i=1:2], SDDP.State, initial_value = i)
    @constraint(subproblem, x[1].in == x[2].out)
    @stageobjective(subproblem, 2.0 * x[1].out)
end

write_to_file(model, joinpath(@__DIR__, "simple.sof.json"))
