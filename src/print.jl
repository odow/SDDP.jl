#  Copyright 2017-20, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function print_helper(f, io, args...)
    f(stdout, args...)
    f(io, args...)
end

function print_banner(io)
    println(
        io,
        "--------------------------------------------------------------------------------",
    )
    println(io, "                      SDDP.jl (c) Oscar Dowson, 2017-20")
    println(io)
end

function print_iteration_header(io)
    println(
        io,
        " Iteration    Simulation       Bound         Time (s)    Proc. ID   # Solves",
    )
end

print_value(x::Real) = lpad(Printf.@sprintf("%1.6e", x), 13)
print_value(x::Int) = Printf.@sprintf("%9d", x)

function print_iteration(io, log::Log)
    print(io, print_value(log.iteration))
    print(io, "   ", print_value(log.simulation_value))
    print(io, "  ", print_value(log.bound))
    print(io, "  ", print_value(log.time))
    print(io, "  ", print_value(log.pid))
    print(io, "  ", print_value(log.total_solves))
    println(io)
end

function print_footer(io, training_results)
    println(io, "\nTerminating training with status: $(training_results.status)")
    println(
        io,
        "------------------------------------------------------------------------------",
    )
end

###
### Numerical stability checks
###

struct CoefficientRanges
    matrix::Vector{Float64}
    objective::Vector{Float64}
    bounds::Vector{Float64}
    rhs::Vector{Float64}
    CoefficientRanges() = new([Inf, -Inf], [Inf, -Inf], [Inf, -Inf], [Inf, -Inf])
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
end

function _update_range(range::Vector{Float64}, func::MOI.LessThan)
    _update_range(range, func.upper)
end

function _update_range(range::Vector{Float64}, func::MOI.GreaterThan)
    _update_range(range, func.lower)
end

function _update_range(range::Vector{Float64}, func::MOI.EqualTo)
    _update_range(range, func.value)
end

function _update_range(range::Vector{Float64}, func::MOI.Interval)
    _update_range(range, func.upper)
    _update_range(range, func.lower)
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
    numerical_stability_report([io::IO=stdout,] model::PolicyGraph,
                               by_node::Bool=false, print=true, warn::Bool=true)

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
    numerical_stability_report(stdout, model, by_node = by_node, print = print, warn = warn)
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
        error("Unable to write the log to file because the model has not been trained.")
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
end
