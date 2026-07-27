# Copyright (c) 2017-26: Oscar Dowson and SDDP.jl contributors.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

_randstring() = joinpath(tempdir(), string(Random.randstring(), ".html"))

"""
    plot(model::PolicyGraph[, filename::String]; open::Bool = true)

This is an experimental function that plots the graph using Javascript.
"""
function plot(
    model::PolicyGraph,
    filename::String = _randstring();
    open::Bool = true,
)
    data = Any[]
    push!(data, "{data: {id: '$(model.root_node)', shape: 'ellipse'}}")
    names = sort(collect(keys(model.nodes)))
    for name in names
        n_terms = length(model[name].noise_terms)
        meta = "Node: $(name)\\nNoise terms: $(n_terms)"
        if n_terms > 1
            push!(data, "{data: {id: '$(name)', meta: '$meta', has_noise: true}}")
        else
            push!(data, "{data: {id: '$(name)', meta: '$meta'}}")
        end
    end
    edge_id = 0
    for child in model.root_children
        edge_id += 1
        meta = "From: $(model.root_node)\\nTo: $(child.term)\\nProbablity: $(child.probability)"
        push!(
            data,
            "{data: {id: 'edge_$(edge_id)', source: '$(model.root_node)', target: '$(child.term)', meta: '$meta'}}",
        )
    end
    for name in names
        for child in model[name].children
            edge_id += 1
            meta = "From: $(name)\\nTo: $(child.term)\\nProbablity: $(child.probability)"
            push!(
                data,
                "{data: {id: 'edge_$(edge_id)', source: '$name', target: '$(child.term)', meta: '$meta'}}",
            )
        end
    end
    fill_template(
        filename,
        "<!--DATA-->" => join(data, ",\n");
        template = joinpath(dirname(@__FILE__), "graph.html"),
        launch = open,
    )
    return
end
