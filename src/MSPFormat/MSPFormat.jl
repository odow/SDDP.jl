#  Copyright (c) 2023, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module MSPFormat

import JSON
import JuMP
import SDDP
import Test

const Object = Dict{String,Any}

_func_variable(name) = Object("type" => "Variable", "name" => name)

function _to_set(type, value)
    if type == "EQ"
        return Object("type" => "EqualTo", "value" => value)
    elseif type == "LEQ"
        return Object("type" => "LessThan", "upper" => value)
    else
        @assert type == "GEQ"
        return Object("type" => "GreaterThan", "lower" => value)
    end
end

function _parse_lattice(filename::String)
    data = JuMP.MOI.FileFormats.compressed_open(
        filename,
        "r",
        JuMP.MOI.FileFormats.AutomaticCompression(),
    ) do io
        return JSON.parse(io)
    end
    graph = SDDP.Graph("root")
    for key in keys(data)
        SDDP.add_node(graph, key)
    end
    for (key, value) in data
        for (child, probability) in value["successors"]
            SDDP.add_edge(graph, key => child, probability)
        end
    end
    stage_zero = String[key for (key, value) in data if value["stage"] == 0]
    for key in stage_zero
        SDDP.add_edge(graph, "root" => key, 1 / length(stage_zero))
    end
    return graph, data
end

_get_constant(terms::String, state::Dict) = get(state, terms, 0.0)

function _get_constant(terms::Vector, state::Union{Dict,Nothing} = nothing)
    if length(terms) == 1
        if terms[1] isa Number
            return terms[1]
        elseif terms[1] == "inf" || terms[1] == "-inf"
            return nothing
        end
    end
    constant = 0.0
    multipliers = 1.0
    for term in terms
        @assert term isa Dict
        key = get(term, "ADD", 0.0)
        if state == nothing && key isa String
            return terms
        end
        value = get(state, key, key)
        if !(value isa Number)
            value = _get_constant(value, state)
        end
        constant += value
        key = get(term, "MUL", 1.0)
        if state == nothing && key isa String
            return terms
        end
        value = get(state, key, key)
        if !(value isa Number)
            value = _get_constant(value, state)
        end
        constant *= value
    end
    return multipliers * constant
end

function test_get_constant()
    Test.@test _get_constant(Any[1.0]) == 1.0
    Test.@test _get_constant(Any[2.4]) == 2.4
    Test.@test _get_constant(Any["inf"]) == nothing
    Test.@test _get_constant(Any["-inf"]) == nothing
    state = Dict{String,Any}("inflow" => 12.0)
    Test.@test _get_constant(Any[1.0], state) == 1.0
    Test.@test _get_constant(Any[2.4], state) == 2.4
    Test.@test _get_constant(Any["inf"], state) == nothing
    Test.@test _get_constant(Any["-inf"], state) == nothing
    Test.@test _get_constant(Any["inf"]) == nothing
    terms = [Dict("ADD" => "inflow"), Dict("ADD" => 0.0)]
    Test.@test _get_constant(terms) === terms
    Test.@test _get_constant(terms, state) === 12.0
    terms = [Dict("ADD" => "inflow"), Dict("ADD" => 200.0), Dict("ADD" => Any[0.0])]
    Test.@test _get_constant(terms) === terms
    Test.@test _get_constant(terms, state) === 212.0
    terms = [Dict("ADD" => "inflow"), Dict("ADD" => 200.0), Dict("ADD" => Any[1.0])]
    Test.@test _get_constant(terms) === terms
    Test.@test _get_constant(terms, state) === 213.0
    return
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

_set_type(::Any, type) = _set_type(0.0, type)

function _build_lhs(stage::Int, sp::JuMP.Model, terms::Vector{Any})
    if maximum(term["stage"] for term in terms) != stage
        return nothing, nothing
    end
    lhs = JuMP.AffExpr(0.0)
    lhs_data = Dict{JuMP.VariableRef,Any}()
    for term in terms
        x = sp[Symbol(term["name"])]
        if x isa JuMP.VariableRef
            @assert term["stage"] == stage
        else
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
            lhs_data[x] = coef
            lhs += x
        else
            lhs += coef * x
        end
    end
    return lhs, lhs_data
end

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

function read_from_file(problem_name::String)
    problem_filename = problem_name * ".problem.json"
    lattice_filename = problem_name * ".lattice.json"
    if !isfile(lattice_filename)
        lattice_filename *= ".gz"
    end
    return read_from_file(problem_filename, lattice_filename)
end

function read_from_file(problem_filename::String, lattice_filename::String)
    graph, graph_data = _parse_lattice(lattice_filename)
    problem = JSON.parsefile(problem_filename)
    state_variables = _state_variables(problem)
    model = SDDP.PolicyGraph(
        graph;
        sense = problem["maximize"] ? :Max : :Min,
        lower_bound = problem["maximize"] ? -Inf : 0,
        upper_bound = problem["maximize"] ? 1e7 : Inf,
    ) do sp, node
        ω_lower_bound = Dict{JuMP.VariableRef,Any}()
        ω_upper_bound = Dict{JuMP.VariableRef,Any}()
        ω_objective = Dict{JuMP.VariableRef,Any}()
        ω_lhs_coefficient = Dict{JuMP.ConstraintRef,Any}()
        ω_rhs_coefficient = Dict{JuMP.ConstraintRef,Any}()
        stage = graph_data[node]["stage"]
        stage_objective = 0.0
        for variable in problem["variables"]
            if variable["stage"] != stage
                continue
            end
            @assert variable["type"] == "CONTINUOUS"
            lower_bound = _get_constant(variable["lb"])
            upper_bound = _get_constant(variable["ub"])
            objective = _get_constant(variable["obj"])
            sym_name = Symbol(variable["name"])
            x = if variable["name"] in state_variables
                x = sp[sym_name] = JuMP.@variable(
                    sp,
                    variable_type = SDDP.State,
                    # Because MSPFormat ignores initial values
                    initial_value = 0.0,
                    base_name = "$sym_name",
                )
                x.out
            else
                sp[sym_name] = JuMP.@variable(sp, base_name = "$sym_name")
            end
            if lower_bound isa Number
                JuMP.set_lower_bound(x, lower_bound)
            elseif lower_bound isa Vector{Any}
                ω_lower_bound[x] = lower_bound
            end
            if upper_bound isa Number
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
        if isempty(stage_objective)
            SDDP.@stageobjective(sp, stage_objective)
        end
        SDDP.parameterize(sp, [graph_data[node]["state"]]) do ω
            SDDP.@stageobjective(
                sp,
                stage_objective +
                sum(x * _get_constant(terms, ω) for (x, terms) in ω_objective)
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

function convert_to_stochoptformat(
    problem_filename::String,
    lattice_filename::String,
    output_filename::Union{Nothing,String} = nothing,
)
    lattice_data = JuMP.MOI.FileFormats.compressed_open(
        lattice_filename,
        "r",
        JuMP.MOI.FileFormats.AutomaticCompression(),
    ) do io
        return JSON.parse(io)
    end
    problem_data = JuMP.MOI.FileFormats.compressed_open(
        problem_filename,
        "r",
        JuMP.MOI.FileFormats.AutomaticCompression(),
    ) do io
        return JSON.parse(io)
    end
    state_variables = _state_variables(problem_data)
    model = Object(
        "version" => Object("major" => 0, "minor" => 2),
        "root" => Object(
            "state_variables" => Object(
                x => Object("initial_value" => 0.0) for
                x in state_variables
            ),
            "successors" => Object(),
        ),
        "nodes" => Object(),
        "subproblems" => Object(),
    )
    for (key, value) in lattice_data
        model["nodes"][key] = Dict(
            "subproblem" => string(value["stage"]),
            "realizations" => [
                Dict("probability" => 1.0, "support" => value["state"])
            ],
            "successors" => value["successors"],
        )
        if value["stage"] == 0
            model["root"]["successors"][key] = 1.0
        end
    end
    for key in keys(model["root"]["successors"])
        model["root"]["successors"][key] /= length(model["root"]["successors"])
    end
    for node in values(model["nodes"])
        if haskey(model["subproblems"], node["subproblem"])
            continue
        end
        model["subproblems"][node["subproblem"]] = Object(
            "state_variables" => Dict(
                x => Dict("in" => "$(x)_in", "out" => x) for
                x in state_variables
            ),
            "random_variables" => String[],
            "subproblem" => Object(
                "version" => Object("major" => 1, "minor" => 2),
                "variables" => Object[],
                "objective" => Object(
                    "sense" => problem_data["maximize"] ? "max" : "min",
                    "function" => Object(
                        "type" => "ScalarQuadraticFunction",
                        "affine_terms" => Object[],
                        "quadratic_terms" => Object[],
                        "constant" => 0.0,
                    )
                ),
                "constraints" => Object[],
            )
        )
    end
    for variable in problem_data["variables"]
        sp = model["subproblems"][string(variable["stage"])]["subproblem"]
        push!(sp["variables"], Object("name" => variable["name"]))
        if variable["name"] in state_variables
            push!(
                sp["variables"],
                Object("name" => variable["name"] * "_in"),
            )
        end
        lb = _get_constant(variable["lb"])
        if lb isa Number
            push!(
                sp["constraints"],
                Object(
                    "function" => _func_variable(variable["name"]),
                    "set" => Object("type" => "GreaterThan", "lower" => lb),
                ),
            )
        elseif lb !== nothing
            @assert false "lower_bound"
        end
        ub = _get_constant(variable["ub"])
        if ub isa Number
            push!(
                sp["constraints"],
                Object(
                    "function" => _func_variable(variable["name"]),
                    "set" => Object("type" => "LessThan", "upper" => ub),
                ),
            )
        elseif ub !== nothing
            @assert false "upper_bound"
        end
        obj = _get_constant(variable["obj"])
        if obj isa Number
            push!(
                sp["objective"]["function"]["affine_terms"],
                Object(
                    "variable" => variable["name"],
                    "coefficient" => obj,
                ),
            )
        elseif obj !== nothing
            @assert false "objective"
        end
    end
    # Restore objective function to affine if not stochastic
    for data in values(model["subproblems"])
        func = data["subproblem"]["objective"]["function"]
        if isempty(func["quadratic_terms"])
            func["type"] = "ScalarAffineFunction"
            func["terms"] = func["affine_terms"]
            delete!(func, "affine_terms")
            delete!(func, "quadratic_terms")
        end
    end
    for constraint in problem_data["constraints"]
        stage = maximum(term["stage"] for term in constraint["lhs"])
        subproblem = model["subproblems"]["$stage"]
        sp = subproblem["subproblem"]
        func = Object(
            "type" => "ScalarQuadraticFunction",
            "affine_terms" => Object[],
            "quadratic_terms" => Object[],
            "constant" => 0.0
        )
        for term in constraint["lhs"]
            x = term["name"]
            if term["stage"] != stage  # Must be an incoming state variable
                x *= "_in"
            end
            coefficient = _get_constant(term["coefficient"])
            if coefficient isa Number
                push!(
                    func["affine_terms"],
                    Object("variable" => x, "coefficient" => coefficient),
                )
            else
                @assert false "lhs"
            end
        end
        rhs = _get_constant(constraint["rhs"])
        set = if rhs isa Number
            _to_set(constraint["type"], rhs)
        else
            rhs_const = 0.0
            for term in rhs
                @assert term isa Dict
                add = term["ADD"]
                if !(add isa String || add isa Number)
                    # Some nested thing?
                    add = _get_constant(add)
                end
                if add isa String
                    if !(add in subproblem["random_variables"])
                        push!(subproblem["random_variables"], add)
                        push!(sp["variables"], Object("name" => add))
                    end
                    push!(
                        func["affine_terms"],
                        Object("variable" => add, "coefficient" => -1.0),
                    )
                elseif !iszero(add)
                    rhs_const += add
                end
            end
            _to_set(constraint["type"], rhs_const)
        end
        # Restore constraint function to affine if not stochastic
        if isempty(func["quadratic_terms"])
            func["type"] = "ScalarAffineFunction"
            func["terms"] = func["affine_terms"]
            delete!(func, "affine_terms")
            delete!(func, "quadratic_terms")
        end
        push!(sp["constraints"], Object("function" => func, "set" => set))
    end
    if output_filename === nothing
        return model
    end
    open(output_filename, "w") do io
        write(io, JSON.json(model))
    end
    return
end

end  # module

