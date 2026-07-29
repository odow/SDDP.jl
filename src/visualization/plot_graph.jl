# Copyright (c) 2017-26: Oscar Dowson and SDDP.jl contributors.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

"""
    plot(model::PolicyGraph[, filename::String]; open::Bool = true)
    plot(model::Graph[, filename::String]; open::Bool = true)

This is an experimental function that plots the structure of the policy graph
`model` using Javascript.

If `filename` is not given, it will be saved to a temporary directory.

If `open = true`, then a browser window will be opened to display the resulting
HTML file.

## How the layout is decided.

This function uses Cytoscape.js to layout the graph. Graph layout is a
notoriously tricky problem. If you have an example where the layout is poor,
please open a GitHub issue and we can experiment with different options in
Cytoscape.
"""
function plot(
    model::PolicyGraph,
    filename::String = joinpath(tempdir(), Random.randstring() * ".html");
    open::Bool = true,
)
    data = Any[]
    meta = join(("$k = $v" for (k, v) in model.initial_root_state), "\\n")
    push!(
        data,
        "{data: {id: '$(model.root_node)', shape: 'ellipse', meta: '$meta'}}",
    )
    names = sort(collect(keys(model.nodes)))
    for name in names
        n_terms = length(model[name].noise_terms)
        meta = "Node: $(name)\\nNoise terms: $(n_terms)"
        if n_terms > 1
            push!(
                data,
                "{data: {id: '$(name)', meta: '$meta', has_noise: true}}",
            )
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

function plot(
    graph::Graph,
    filename::String = joinpath(tempdir(), Random.randstring() * ".html");
    open::Bool = true,
)
    data = Any[]
    push!(data, "{data: {id: '$(graph.root_node)', shape: 'ellipse'}}")
    names = sort(collect(keys(model.nodes)))
    for name in names
        push!(data, "{data: {id: '$(name)', meta: 'Node: $name'}}")
    end
    edge_id = 0
    for name in vcat(graph.root_node, names)
        for (child, probability) in model[name]
            edge_id += 1
            meta = "From: $(name)\\nTo: $child\\nProbablity: $probability"
            push!(
                data,
                "{data: {id: 'edge_$(edge_id)', source: '$name', target: '$child', meta: '$meta'}}",
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
