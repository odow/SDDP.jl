#  Copyright (c) 2017-22, Oscar Dowson and SDDP.jl contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function print_helper(f, io, args...)
    f(stdout, args...)
    return f(io, args...)
end

function print_banner(io)
    println(
        io,
        "------------------------------------------------------------------------------",
    )
    println(io, "          SDDP.jl (c) Oscar Dowson and SDDP.jl contributors, 2017-21")
    return println(io)
end

function _unique_paths(model::PolicyGraph{T}) where {T}
    if is_cyclic(model)
        return Inf
    end
    parents = Dict{T,Set{T}}(t => Set{T}() for t in keys(model.nodes))
    children = Dict{T,Set{T}}(t => Set{T}() for t in keys(model.nodes))
    for (t, node) in model.nodes
        for child in node.children
            if child.probability > 0
                push!(parents[child.term], t)
                push!(children[t], child.term)
            end
        end
    end
    ordered = T[]
    in_order = Dict{T,Bool}(t => false for t in keys(model.nodes))
    stack = Tuple{T,Bool}[]
    for root_child in model.root_children
        if iszero(root_child.probability) || in_order[root_child.term]
            continue
        end
        push!(stack, (root_child.term, true))
        while !isempty(stack)
            node, needs_checking = pop!(stack)
            if !needs_checking
                push!(ordered, node)
                in_order[node] = true
                continue
            elseif in_order[node]
                continue
            end
            push!(stack, (node, false))
            for child in children[node]
                if !in_order[child]
                    push!(stack, (child, true))
                end
            end
        end
    end
    total_scenarios = 0.0
    incoming_scenarios = Dict{T,Float64}(t => 0.0 for t in keys(model.nodes))
    for node in reverse!(ordered)
        N = length(model[node].noise_terms)
        if length(parents[node]) == 0  # Must come from the root node.
            incoming_scenarios[node] = N
        else
            incoming_scenarios[node] =
                N * sum(incoming_scenarios[p] for p in parents[node])
        end
        if length(children[node]) == 0  # It's a leaf!
            total_scenarios += incoming_scenarios[node]
        end
    end
    return total_scenarios
end

function _merge_tuple(x, y)
    if x == (-1, -1)
        return (y, y)
    elseif y < x[1]
        return (y, x[2])
    elseif y > x[2]
        return (x[1], y)
    else
        return x
    end
end

_constraint_key(F, S) = replace("$(F) in $(S)", "MathOptInterface" => "MOI")

function print_problem_statistics(
    io::IO,
    model::PolicyGraph,
    existing_cuts::Bool,
    parallel_scheme,
    risk_measure,
    sampling_scheme,
)
    constraint_types = Dict{String,Tuple{Int,Int}}()
    variables = (-1, -1)
    for (_, node) in model.nodes
        variables = _merge_tuple(variables, JuMP.num_variables(node.subproblem))
        for (F, S) in JuMP.list_of_constraint_types(node.subproblem)
            key = _constraint_key(F, S)
            num_con = get(constraint_types, key, (-1, -1))
            constraint_types[key] = _merge_tuple(
                num_con,
                JuMP.num_constraints(node.subproblem, F, S),
            )
        end
    end
    pad = maximum(length(k) for k in keys(constraint_types))
    println(io, "Problem")
    println(io, "  Nodes           : ", length(model.nodes))
    println(io, "  State variables : ", length(model.initial_root_state))
    paths = Printf.@sprintf("%1.5e", _unique_paths(model))
    println(io, "  Scenarios       : ", paths)
    println(io, "  Existing cuts   : ", existing_cuts)
    println(io, rpad("  Subproblem structure", pad + 4), " : (min, max)")
    println(io, "    ", rpad("Variables", pad), " : ", variables)
    for (k, v) in constraint_types
        println(io, "    ", rpad(k, pad), " : ", v)
    end
    println(io, "Options")
    println(io, "  Solver          : ", parallel_scheme)
    println(io, "  Risk measure    : ", risk_measure)
    println(io, "  Sampling scheme : ", typeof(sampling_scheme))
    println(io)
    return
end

function print_iteration_header(io)
    println(
        io,
        " Iteration    Simulation       Bound         Time (s)    Proc. ID   # Solves",
    )
    return
end

print_value(x::Real) = lpad(Printf.@sprintf("%1.6e", x), 13)
print_value(x::Int) = Printf.@sprintf("%9d", x)

function print_iteration(io, log::Log)
    print(io, print_value(log.iteration))
    print(io, log.duality_key)
    print(io, "  ", print_value(log.simulation_value))
    print(io, "  ", print_value(log.bound))
    print(io, "  ", print_value(log.time))
    print(io, "  ", print_value(log.pid))
    print(io, "  ", print_value(log.total_solves))
    println(io)
    return
end

function print_footer(io, training_results::TrainingResults)
    println(io)
    println(io, "Terminating training")
    println(io, "  Status         : ", training_results.status)
    println(
        io,
        "  Total time (s) :",
        print_value(training_results.log[end].time),
    )
    println(io, "  Total solves   : ", training_results.log[end].total_solves)
    println(
        io,
        "  Best bound     : ",
        print_value(training_results.log[end].bound),
    )
    μ, σ =
        confidence_interval(map(l -> l.simulation_value, training_results.log))
    println(io, "  Simulation CI  : ", print_value(μ), " ±", print_value(σ))
    println(
        io,
        "------------------------------------------------------------------------------",
    )
    return
end

"""
    confidence_interval(x::Vector{Float64}, z_score::Float64 = 1.96)

Return a confidence interval of `x` corresponding to the `z_score`.

`z_score` defaults to `1.96` for a 95% confidence interval.
"""
function confidence_interval(x::Vector{Float64}, z_score::Float64 = 1.96)
    μ = Statistics.mean(x)
    σ = z_score * Statistics.std(x) / sqrt(length(x))
    return μ, σ
end

###
### Numerical stability checks
###

struct CoefficientRanges
    matrix::Vector{Float64}
    objective::Vector{Float64}
    bounds::Vector{Float64}
    rhs::Vector{Float64}
    function CoefficientRanges()
        return new([Inf, -Inf], [Inf, -Inf], [Inf, -Inf], [Inf, -Inf])
    end
end

function _merge(x::Vector{Float64}, y::Vector{Float64})
    x[1] = min(x[1], y[1])
    x[2] = max(x[2], y[2])
    return
end
function _merge(x::CoefficientRanges, y::CoefficientRanges)
    _merge(x.matrix, y.matrix)
    _merge(x.objective, y.objective)
    _merge(x.bounds, y.bounds)
    _merge(x.rhs, y.rhs)
    return
end

function _stringify_bounds(bounds::Vector{Float64})
    lower = bounds[1] < Inf ? _print_value(bounds[1]) : "0e+00"
    upper = bounds[2] > -Inf ? _print_value(bounds[2]) : "0e+00"
    return string("[", lower, ", ", upper, "]")
end

function _print_numerical_stability_report(
    io::IO,
    ranges::CoefficientRanges,
    print::Bool,
    warn::Bool,
)
    warnings = Tuple{String,String}[]
    _print_coefficients(io, "Matrix", ranges.matrix, print, warnings)
    _print_coefficients(io, "Objective", ranges.objective, print, warnings)
    _print_coefficients(io, "Bounds", ranges.bounds, print, warnings)
    _print_coefficients(io, "RHS", ranges.rhs, print, warnings)
    if warn
        if length(warnings) > 0
            println(io, "WARNING: numerical stability issues detected")
            for (name, sense) in warnings
                println(io, "  - $(name) range contains $(sense) coefficients")
            end
            println(
                io,
                "Very large or small absolute values of coefficients\n",
                "can cause numerical stability issues. Consider\n",
                "reformulating the model.",
            )
        else
            print && println(io, "No problems detected")
        end
    end
    return
end

function _print_coefficients(
    io::IO,
    name::String,
    range,
    print::Bool,
    warnings::Vector{Tuple{String,String}},
)
    if print
        println(
            io,
            "  Non-zero ",
            rpad(string(name, " range"), 17),
            _stringify_bounds(range),
        )
    end
    range[1] < 1e-4 && push!(warnings, (name, "small"))
    range[2] > 1e7 && push!(warnings, (name, "large"))
    return
end

_print_value(x::Real) = Printf.@sprintf("%1.0e", x)

function _update_range(range::Vector{Float64}, value::Real)
    if !(value ≈ 0.0)
        range[1] = min(range[1], abs(value))
        range[2] = max(range[2], abs(value))
    end
    return
end

function _update_range(range::Vector{Float64}, func::JuMP.GenericAffExpr)
    for coefficient in values(func.terms)
        _update_range(range, coefficient)
    end
    return
end

function _update_range(range::Vector{Float64}, func::MOI.LessThan)
    _update_range(range, func.upper)
    return
end

function _update_range(range::Vector{Float64}, func::MOI.GreaterThan)
    _update_range(range, func.lower)
    return
end

function _update_range(range::Vector{Float64}, func::MOI.EqualTo)
    _update_range(range, func.value)
    return
end

function _update_range(range::Vector{Float64}, func::MOI.Interval)
    _update_range(range, func.upper)
    _update_range(range, func.lower)
    return
end

# Default fallback for unsupported constraints.
_update_range(range::Vector{Float64}, x) = nothing

function _coefficient_ranges(model::JuMP.Model)
    ranges = CoefficientRanges()
    _update_range(ranges.objective, JuMP.objective_function(model))
    for var in JuMP.all_variables(model)
        if JuMP.has_lower_bound(var)
            _update_range(ranges.bounds, JuMP.lower_bound(var))
        end
        if JuMP.has_upper_bound(var)
            _update_range(ranges.bounds, JuMP.upper_bound(var))
        end
    end
    for (F, S) in JuMP.list_of_constraint_types(model)
        F == JuMP.VariableRef && continue
        for con in JuMP.all_constraints(model, F, S)
            con_obj = JuMP.constraint_object(con)
            _update_range(ranges.matrix, con_obj.func)
            _update_range(ranges.rhs, con_obj.set)
        end
    end
    return ranges
end

"""
    numerical_stability_report(
        [io::IO=stdout,]
        model::PolicyGraph,
        by_node::Bool = false,
        print::Bool = true,
        warn::Bool = true,
    )

Print a report identifying possible numeric stability issues.

- If `by_node`, print a report for each node in the graph.
- If `print`, print to `io`.
- If `warn`, warn if the coefficients may cause numerical issues.
"""
function numerical_stability_report(
    io::IO,
    model::PolicyGraph;
    by_node::Bool = false,
    print::Bool = true,
    warn::Bool = true,
)
    graph_ranges = CoefficientRanges()
    node_keys = sort_nodes(collect(keys(model.nodes)))
    for key in node_keys
        node = model[key]
        node_ranges = CoefficientRanges()
        for noise in node.noise_terms
            parameterize(node, noise.term)
            node_ranges_2 = _coefficient_ranges(node.subproblem)
            _merge(node_ranges, node_ranges_2)
        end
        if by_node
            print && println(io, "Numerical stability report for node: ", key)
            _print_numerical_stability_report(io, node_ranges, print, warn)
        end
        _merge(graph_ranges, node_ranges)
    end
    if !by_node
        print && println(io, "Numerical stability report")
        _print_numerical_stability_report(io, graph_ranges, print, warn)
    end
    print && println(io)
    return
end

function numerical_stability_report(
    model::PolicyGraph;
    by_node::Bool = false,
    print::Bool = true,
    warn::Bool = true,
)
    return numerical_stability_report(
        stdout,
        model,
        by_node = by_node,
        print = print,
        warn = warn,
    )
end

###
### Machine readable log
###

"""
    write_log_to_csv(model::PolicyGraph, filename::String)

Write the log of the most recent training to a csv for post-analysis.

Assumes that the model has been trained via [`SDDP.train`](@ref).
"""
function write_log_to_csv(model::PolicyGraph, filename::String)
    if model.most_recent_training_results === nothing
        error(
            "Unable to write the log to file because the model has not " *
            "been trained.",
        )
    end
    open(filename, "w") do io
        println(io, "iteration, simulation, bound, time")
        for log in model.most_recent_training_results.log
            println(
                io,
                log.iteration,
                ", ",
                log.simulation_value,
                ", ",
                log.bound,
                ", ",
                log.time,
            )
        end
    end
    return
end
