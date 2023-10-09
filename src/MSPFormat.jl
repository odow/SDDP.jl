#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module MSPFormat

import JSON
import JuMP
import ..SDDP

function _parse_lattice(filename::String)
    data = JuMP.MOI.FileFormats.compressed_open(
        JSON.parse,
        filename,
        "r",
        JuMP.MOI.FileFormats.AutomaticCompression(),
    )
    graph = SDDP.Graph("root")
    for key in keys(data)
        SDDP.add_node(graph, key)
    end
    for (key, value) in data
        for (child, probability) in sort(value["successors"])
            SDDP.add_edge(graph, key => child, probability)
        end
    end
    # MSPFormat doesn't have explicit root -> stage 1 arcs. Assume uniform.
    # Also, MSPFormat uses 0-indexed stages.
    stage_zero = String[key for (key, value) in data if value["stage"] == 0]
    for key in stage_zero
        SDDP.add_edge(graph, "root" => key, 1 / length(stage_zero))
    end
    return _reduce_lattice(graph, data)
end

"""
    _reduce_lattice(graph, data)

This function takes a graph and associated data from `_parse_lattice`, which is
assumed to represent a Markovian lattice, and tries to collapse nodes in the
graph to be stagewise independent (if possible).

For now, a trivial way to detect stagewise independence is to group nodes by
their unique `graph.nodes`, and then check that the nodes in the group match the
nodes in a stage.
"""
function _reduce_lattice(graph, data)
    arcs_to_node = Dict{UInt64,Vector{String}}()
    for (node, arcs) in graph.nodes
        push!(get!(() -> String[], arcs_to_node, hash(arcs)), node)
    end
    for v in values(arcs_to_node)
        sort!(v)
    end
    nodes_by_stage = Dict{Int,Vector{String}}()
    for (node, d) in data
        push!(get!(() -> String[], nodes_by_stage, d["stage"]), node)
    end
    for v in values(nodes_by_stage)
        sort!(v)
    end
    if !all(n -> n in Set(values(arcs_to_node)), values(nodes_by_stage))
        # Model is not stagewise independent
        graph_data = Dict(
            k => Dict(
                "stage" => v["stage"],
                "sample_space" => [v["state"]],
                "probability" => [1.0],
            ) for (k, v) in data
        )
        return graph, graph_data
    end
    new_graph = SDDP.Graph("root")
    for t in keys(nodes_by_stage)
        SDDP.add_node(new_graph, "$t")
    end
    for t in keys(nodes_by_stage)
        if t == 0
            SDDP.add_edge(new_graph, "root" => "$t", 1.0)
        else
            SDDP.add_edge(new_graph, "$(t-1)" => "$t", 1.0)
        end
    end
    graph_data = Dict(
        "$t" => Dict{String,Any}(
            "stage" => t,
            "sample_space" => Any[],
            "probability" => Float64[],
        ) for t in keys(nodes_by_stage)
    )
    parent = "root"
    while !isempty(graph.nodes[parent])
        for (node, probability) in sort(graph.nodes[parent])
            t = string(data[node]["stage"])
            push!(graph_data[t]["sample_space"], data[node]["state"])
            push!(graph_data[t]["probability"], probability)
            parent = node
        end
    end
    return new_graph, graph_data
end

# Use a default of 0.0 for any missing keys.
_get_constant(terms::String, state::Dict) = get(state, terms, 0.0)
_get_constant(::String, ::Nothing) = nothing
_get_constant(key::Number, ::Union{Dict,Nothing}) = key

function _get_constant(terms::Vector, state::Union{Dict,Nothing} = nothing)
    if length(terms) == 1
        # Special case: if `terms = Any[1.0]` or `terms = Any["inf"]`, then we
        # don't need complicated recursive logic. Bail early.
        if terms[1] isa Number
            return terms[1]
        elseif terms[1] == "inf"
            return Inf
        elseif terms[1] == "-inf"
            return -Inf
        elseif terms[1] isa String
            value = _get_constant(terms[1], state)
            return something(value, terms)
        end
    end
    result = nothing
    for term in terms
        @assert term isa Dict
        if haskey(term, "ADD")
            value = _get_constant(term["ADD"], state)
            if value === nothing
                return terms
            end
            result = something(result, 0.0) + value
        else
            @assert haskey(term, "MUL")
            value = _get_constant(term["MUL"], state)
            if value === nothing
                return terms
            end
            result = something(result, 1.0) * value
        end
    end
    return result::Number
end

function _set_type(rhs::Number, type)
    if type == "EQ"
        return JuMP.MOI.EqualTo{Float64}(rhs)
    elseif type == "LEQ"
        return JuMP.MOI.LessThan{Float64}(rhs)
    else
        @assert type == "GEQ"
        return JuMP.MOI.GreaterThan{Float64}(rhs)
    end
end

# If the RHS is not a Number, it must be a random expression. Use a default RHS
# of 0.0.
_set_type(::Any, type) = _set_type(0.0, type)

function _build_lhs(stage::Integer, sp::JuMP.Model, terms::Vector{Any})
    if maximum(term["stage"] for term in terms) != stage
        # Skip constraints which are not relevant for this stage.
        return nothing, nothing
    end
    # For now, we assume the LHS is affine.
    lhs = JuMP.AffExpr(0.0)
    # lhs_data will store random coefficient terms for each variable.
    lhs_data = Dict{JuMP.VariableRef,Any}()
    for term in terms
        # Lookup variable by name from the JuMP model.
        x = sp[Symbol(term["name"])]
        if x isa JuMP.VariableRef
            @assert term["stage"] == stage
        else
            # `x` is a state, so we need to distinguish whether we want the
            # `.in` or `.out` variables.
            @assert x isa SDDP.State
            if term["stage"] == stage
                x = x.out
            else
                @assert term["stage"] == stage - 1
                x = x.in
            end
        end
        coef = _get_constant(term["coefficient"])
        if coef isa Vector{Any}
            lhs_data[x] = coef  # Store the random data
            # Set lhs += 1.0 * x for now. This will get updated in parameterize.
            lhs += x
        else
            lhs += coef * x
        end
    end
    return lhs, lhs_data
end

# MSPFormat does not store an explicit list of state variables. Detect state
# variables by finding two variables in the same constraint with the same name
# and different `stage` values.
function _state_variables(problem)
    states = Set{String}()
    for constraint in problem["constraints"]
        terms = constraint["lhs"]
        stage = maximum(term["stage"] for term in terms)
        for term in terms
            if term["stage"] != stage
                push!(states, term["name"])
            end
        end
    end
    return sort(collect(states))
end

"""
    read_from_file(
        problem_filename::String,
        lattice_filename::String;
        bound::Float64 = 1e6,
    )

Return a [`SDDP.PolicyGraph`](@ref) built from the MSPFormat files
`problem_filename` and `lattice_filename`, which point to the `.problem.json`
and `.lattice.json` files respectively.

!!! warning
    This function is experimental and may change in any future commit.

## Keyword arguments

 * `bound::Float64 = 1e6`. The absolute value of the lower bound (if minimizing)
   or the upper bound (if maximizing).
"""
function read_from_file(
    problem_filename::String,
    lattice_filename::String;
    bound::Float64 = 1e6,
)
    graph, graph_data = _parse_lattice(lattice_filename)
    problem = JSON.parsefile(problem_filename)
    state_variables = _state_variables(problem)
    initial_values = Dict{Symbol,Float64}()
    model = SDDP.PolicyGraph(
        graph;
        sense = problem["maximize"] ? :Max : :Min,
        lower_bound = problem["maximize"] ? -Inf : -bound,
        upper_bound = problem["maximize"] ? bound : Inf,
    ) do sp, node
        ω_lower_bound = Dict{JuMP.VariableRef,Any}()
        ω_upper_bound = Dict{JuMP.VariableRef,Any}()
        ω_objective = Dict{JuMP.VariableRef,Any}()
        ω_lhs_coefficient = Dict{JuMP.ConstraintRef,Any}()
        ω_rhs_coefficient = Dict{JuMP.ConstraintRef,Any}()
        stage = graph_data[node]["stage"]
        stage_objective = JuMP.AffExpr(0.0)
        for variable in problem["variables"]
            if variable["stage"] != stage
                continue
            end
            lower_bound = _get_constant(variable["lb"])
            upper_bound = _get_constant(variable["ub"])
            objective = _get_constant(variable["obj"])
            sym_name = Symbol(variable["name"])
            if !haskey(initial_values, sym_name)
                initial_values[sym_name] = 0.0
            end
            if lower_bound isa Number && isfinite(lower_bound)
                initial_values[sym_name] =
                    max(initial_values[sym_name], lower_bound)
            end
            if upper_bound isa Number && isfinite(upper_bound)
                initial_values[sym_name] =
                    min(initial_values[sym_name], upper_bound)
            end
            x = if variable["name"] in state_variables
                sp[sym_name] = JuMP.@variable(
                    sp,
                    variable_type = SDDP.State,
                    initial_value = initial_values[sym_name],
                    base_name = "$sym_name",
                )
                sp[sym_name].out
            else
                sp[sym_name] = JuMP.@variable(sp, base_name = "$sym_name")
            end
            if variable["type"] == "BINARY"
                set_binary(x)
            elseif variable["type"] == "INTEGER"
                set_integer(x)
            else
                @assert variable["type"] == "CONTINUOUS"
            end
            if lower_bound isa Number && isfinite(lower_bound)
                JuMP.set_lower_bound(x, lower_bound)
            elseif lower_bound isa Vector{Any}
                ω_lower_bound[x] = lower_bound
            end
            if upper_bound isa Number && isfinite(upper_bound)
                JuMP.set_upper_bound(x, upper_bound)
            elseif upper_bound isa Vector{Any}
                ω_upper_bound[x] = upper_bound
            end
            if objective isa Number
                stage_objective += objective * x
            elseif objective isa Vector{Any}
                ω_objective[x] = objective
            end
        end
        for str_name in state_variables
            sym_name = Symbol(str_name)
            if !haskey(JuMP.object_dictionary(sp), sym_name)
                sp[sym_name] = JuMP.@variable(
                    sp,
                    variable_type = SDDP.State,
                    initial_value = 0.0,
                    base_name = "$sym_name",
                )
            end
        end
        for constraint in problem["constraints"]
            lhs, lhs_data = _build_lhs(stage, sp, constraint["lhs"])
            if lhs === nothing
                continue
            end
            rhs = _get_constant(constraint["rhs"])
            set = _set_type(rhs, constraint["type"])
            con = JuMP.@constraint(sp, lhs in set)
            if rhs isa Vector{Any}
                ω_rhs_coefficient[con] = rhs
            end
            if lhs_data !== nothing
                ω_lhs_coefficient[con] = lhs_data
            end
        end
        Ω = graph_data[node]["sample_space"]
        P = graph_data[node]["probability"]
        SDDP.parameterize(sp, Ω, P) do ω
            SDDP.@stageobjective(
                sp,
                stage_objective + sum(
                    x * _get_constant(terms, ω) for (x, terms) in ω_objective
                )
            )
            for (x, terms) in ω_lower_bound
                JuMP.set_lower_bound(x, _get_constant(terms, ω))
            end
            for (x, terms) in ω_upper_bound
                JuMP.set_upper_bound(x, _get_constant(terms, ω))
            end
            for (con, lhs_data) in ω_lhs_coefficient
                for (x, terms) in lhs_data
                    JuMP.set_normalized_coefficient(con, x, _get_constant(terms, ω))
                end
            end
            for (con, terms) in ω_rhs_coefficient
                JuMP.set_normalized_rhs(con, _get_constant(terms, ω))
            end
        end
        return
    end
    return model
end

"""
    read_from_file(problem_name::String; kwargs...)

A utility for reading MSPFormat files that saves writing out both the problem
and lattice filenames if they are in the same location and differ only by the
suffix.

It is equivalent to a call like:

```julia
read_from_file(problem_name * ".problem.json", problem_name * ".lattice.json")
```

In addition, this function searches for compressed `.gz` versions of the lattice
file, since it may be very large.
"""
function read_from_file(problem_name::String; kwargs...)
    problem_filename = problem_name * ".problem.json"
    lattice_filename = problem_name * ".lattice.json"
    if !isfile(lattice_filename)
        lattice_filename *= ".gz"
    end
    return read_from_file(problem_filename, lattice_filename; kwargs...)
end

end  # module
