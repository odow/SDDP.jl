# ========================== Begin MIT Licensed Code ========================= #
# The following a modified version of that found at
#
# Humanize.jl    https://github.com/IainNZ/Humanize.jl
# as it was at commit be2c55008b501e17ed13c0a9aa791d40214385ea
#
# Copyright (c) 2017 Oscar Dowson, Iain Dunning, Julian Gehring
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software. # THE SOFTWARE IS PROVIDED "AS
# IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
# AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function humanize(value::Integer)
    return Printf.@sprintf("% 9d", value)
end

function humanize(value::Real)
    SUFFIX = [" ", "K", "M", "G", "T", "P", "E", "Z", "Y"]
    BASE = 1000.0
    bytes = abs(float(value)) # O.D. abs value
    unit = BASE
    suffix = SUFFIX[1]
    for (i, s) in enumerate(SUFFIX)
        unit = BASE ^ i
        suffix = s
        if bytes < unit
            break
        end
    end
    normalized_value = sign(value) * BASE * bytes / unit
    return Printf.@sprintf("% 8.3f%s", normalized_value, suffix)
end
# =========================== End MIT Licensed Code ========================== #

# ======================== Begin MPL2.0 Licensed Code ======================== #
#  Copyright 2018, Oscar Dowson. This Source Code Form is subject to the terms
#  of the Mozilla Public License, v.2.0. If a copy of the MPL was not
#  distributed with this file, You can obtain one at
#  http://mozilla.org/MPL/2.0/.

function print_banner(io=stdout)
    println(io, "———————————————————————————————————————————————————————————————————————————————")
    println(io, "                        SDDP.jl - © Oscar Dowson, 2017-19.")
    println(io, "———————————————————————————————————————————————————————————————————————————————")
end

function print_iteration_header(io=stdout)
    println(io, " Iteration | Simulation |      Bound |   Time (s)")
    println(io, "———————————————————————————————————————————————————————————————————————————————")
end

function print_iteration(io, log::Log)
    line = string(" ", humanize(log.iteration), " |  ",
        humanize(log.simulation_value), " |  ", humanize(log.bound), " |  ",
        humanize(log.time))
    println(io, line)
end

function print_footer(io, training_results)
    println(io, "———————————————————————————————————————————————————————————————————————————————")
    println(io, " Terminating training with status: $(training_results.status)")
    println(io, "———————————————————————————————————————————————————————————————————————————————")
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
        io::IO, ranges::CoefficientRanges, print::Bool, warn::Bool)
    warnings = Tuple{String, String}[]
    _print_coefficients(io, "Matrix", ranges.matrix, print, warnings)
    _print_coefficients(io, "Objective", ranges.objective, print, warnings)
    _print_coefficients(io, "Bounds", ranges.bounds, print, warnings)
    _print_coefficients(io, "RHS", ranges.rhs, print, warnings)
    if warn && length(warnings) > 0
        println(io, "WARNING: numerical stability issues detected")
        for (name, sense) in warnings
            println(io, "  - $(name) range contains $(sense) coefficients")
        end
        println(io, "These coefficients can cause numerical stability issues. ",
                "Consider reformulating\nthe model.")
    end
    return
end

function _print_coefficients(io::IO, name::String, range, print::Bool, warnings::Vector{Tuple{String, String}})
    if print
        println(io, "  Non-zero ", rpad(string(name, " range"), 17),
                _stringify_bounds(range))
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
        io::IO, model::PolicyGraph;
        by_node::Bool=false, print::Bool=true, warn::Bool=true)
    graph_ranges = CoefficientRanges()
    node_keys = sort_nodes(collect(keys(model.nodes)))
    for key in node_keys
        node = model[key]
        node_ranges = CoefficientRanges()
        for noise in node.noise_terms
            node.parameterize(noise.term)
            set_objective(model, node)
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
    print && println(io, "———————————————————————————————————————————————————————————————————————————————")
    return
end

function numerical_stability_report(
        model::PolicyGraph;
        by_node::Bool=false, print::Bool=true, warn::Bool=true)
    numerical_stability_report(
        stdout, model, by_node=by_node, print=print, warn=warn)
end
