#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

import SHA

"""
    TestScenario{T, S}(probability::Float64, scenario::Vector{Tuple{T, S}})

A single scenario for testing.

See also: [`TestScenarios`](@ref).
"""
struct TestScenario{T, S}
    probability::Float64
    scenario::Vector{Tuple{T, S}}
end

"""
    TestScenarios{T, S}(scenarios::Vector{TestScenario{T, S}})

An [`AbstractSamplingScheme`](@ref) based on a vector of scenarios.

Each scenario is a vector of `Tuple{T, S}` where the first element is the node
to visit and the second element is the realization of the stagewise-independent
noise term. Pass `nothing` if the node is deterministic.
"""
mutable struct TestScenarios{T, S} <: AbstractSamplingScheme
    scenarios::Vector{TestScenario{T, S}}
    last::Int
    SHA256::String

    function TestScenarios(
        scenarios::Vector{TestScenario{T, S}}; SHA256::String = ""
    ) where {T, S}
        return new{T, S}(scenarios, 0, SHA256)
    end
end

function sample_scenario(
    model::PolicyGraph{T}, sampling_scheme::TestScenarios{T, S}; kwargs...
) where {T, S}
    sampling_scheme.last += 1
    if sampling_scheme.last > length(sampling_scheme.scenarios)
        sampling_scheme.last = 1
    end
    return sampling_scheme.scenarios[sampling_scheme.last].scenario, false
end

function _throw_if_belief_states(model::PolicyGraph)
    if length(model.belief_partition) != 0
        error("StochOptFormat does not support belief states.")
    end
end

function _throw_if_objective_states(model::PolicyGraph)
    for (_, node) in model.nodes
        if node.objective_state !== nothing
            error("StochOptFormat does not support objective states.")
        end
    end
end

function _throw_if_exisiting_cuts(model::PolicyGraph)
    for (_, node) in model.nodes
        if length(node.bellman_function.global_theta.cut_oracle.cuts) != 0
            error(
                "StochOptFormat does not support writing after a call to " *
                "`SDDP.train`."
            )
        end
    end
end

function _test_scenarios(model::PolicyGraph, test_scenarios::Int, scenario_map)
    return _test_scenarios(
        model,
        TestScenarios([
            TestScenario(
                1 / test_scenarios,
                sample_scenario(model, InSampleMonteCarlo())[1]
            )
            for _ = 1:test_scenarios
        ]),
        scenario_map,
    )
end
function _test_scenarios(
    ::PolicyGraph, test_scenarios::TestScenarios, scenario_map
)
    return [
        Dict(
            "probability" => scenario.probability,
            "scenario" => [
                Dict(
                    "node" => "$(node)",
                    "support" => scenario_map[node][noise]
                ) for (node, noise) in scenario.scenario
            ]
        )
        for scenario in test_scenarios.scenarios
    ]
end

"""
    Base.write(
        io::IO,
        model::PolicyGraph;
        test_scenarios::Union{Int, TestScenarios} = 1_000,
        kwargs...
    )

Write `model` to `io` in the StochOptFormat file format.

Pass an `Int` to `test_scenarios` (default `1_000`) to specify the number of
test scenarios to generate using the [`InSampleMonteCarlo`](@ref) sampling
scheme. Alternatively, pass a [`TestScenarios`](@ref) object to manually specify
the test scenarios to use.

Any additional `kwargs` passed to `write` will be stored in the top-level of the
resulting StochOptFormat file. Valid arguments include `name`, `author`, `date`,
and `description`.

WARNING: THIS FUNCTION IS EXPERIMENTAL. THINGS MAY CHANGE BETWEEN COMMITS. YOU
SHOULD NOT RELY ON THIS FUNCTIONALITY AS A LONG-TERM FILE FORMAT (YET).

In addition to potential changes to the underlying format, only a subset of
possible modifications are supported. These include:
- `JuMP.fix`
- `JuMP.set_lower_bound`
- `JuMP.set_upper_bound`
- `JuMP.set_normalized_rhs`
- Changes to the constant or affine terms in a stage objective

If your model uses something other than this, this function will silently write
an incorrect formulation of the problem.

## Example

    open("my_model.sof.json", "w") do io
        write(
            io,
            model;
            test_scenarios = 10,
            name = "MyModel",
            author = "@odow",
            date = "2020-07-20",
            description = "Example problem for the SDDP.jl documentation",
        )
    end
"""
function Base.write(
    io::IO,
    model::PolicyGraph{T};
    test_scenarios::Union{Int, TestScenarios{T, S}} = 1_000,
    kwargs...
) where {T, S}
    _throw_if_belief_states(model)
    _throw_if_objective_states(model)
    _throw_if_exisiting_cuts(model)
    edges = Dict{String, Any}[]
    _add_edges(edges, "$(model.root_node)", model.root_children)
    nodes = Dict{String, Any}()
    scenario_map = Dict{T, Any}()
    for (node_name, node) in model.nodes
        _add_edges(edges, "$(node_name)", node.children)
        scenario_map[node_name] = _add_node_to_dict(nodes, node, "$(node_name)")
    end
    sof = Dict{String, Any}(
        "version" => Dict("major" => 0, "minor" => 1),
        "root" => Dict{String, Any}(
            "name" => "$(model.root_node)",
            "state_variables" => Dict{String, Any}(
                "$(k)" => Dict{String, Any}("initial_value" => v)
                for (k, v) in model.initial_root_state
            )
        ),
        "nodes" => nodes,
        "edges" => edges,
        "test_scenarios" => _test_scenarios(model, test_scenarios, scenario_map)
    )
    for (k, v) in kwargs
        sof["$(k)"] = v
    end
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
    random_variables = String[]
    realizations = Dict{String, Any}[
        Dict{String, Any}(
            "probability" => noise.probability,
            "support" => Dict{String, Float64}()
        ) for noise in node.noise_terms
    ]
    undo_reformulation = _reformulate_uncertainty(
        node, realizations, random_variables
    )
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
    undo_reformulation()
    # Return a dictionary which maps the in-sample realizations to the
    # transformed support for writing `test_scenarios` to file.
    return Dict(
        noise.term => realizations[i]["support"]
        for (i, noise) in enumerate(node.noise_terms)
    )
end

"""
    _reformulate_uncertainty(
        node::Node, realizations, random_variables
    )

Convert any lower and upper variable_bound_storage than depend on the uncertainty into linear
constraints with a random variable.

Fixed variables are recorded as random variables, but no transformation is done.

Given an affine stageobjective with stagewise independent uncertainty,
reformulate into a quadratic stage objective by replacing the random
coefficients with random decision variables.

Return a function that undoes the reformulation when called with no arguments.
"""
function _reformulate_uncertainty(
    node::Node, realizations::Vector, random_variables::Vector{String}
)
    # Storage for things that are changing.
    variable_bound_storage = Dict{VariableRef, Any}[]
    changing_variable_lower_bound = Set{VariableRef}()
    changing_variable_upper_bound = Set{VariableRef}()
    changing_variable_fixed_bound = Set{VariableRef}()

    objective_storage = AffExpr[]
    changing_objective_constant = false
    changing_objective_coefficient = Set{VariableRef}()

    constraint_rhs_storage = Dict{ConstraintRef, Float64}[]
    changing_constraint_rhs = Set{ConstraintRef}()

    # Collect terms that are changing
    for noise in node.noise_terms
        node.parameterize(noise.term)

        # Collect changing variable bounds.
        _collect_changing_variable_bounds(
            node,
            variable_bound_storage,
            changing_variable_lower_bound,
            changing_variable_upper_bound,
            changing_variable_fixed_bound,
        )

        # Collect changing objective terms.
        changing_objective_constant = _collect_changing_objective(
            node,
            objective_storage,
            changing_objective_constant,
            changing_objective_coefficient,
        )

        # Collect changing RHS terms.
        _collect_changing_constraint_rhs(
            node,
            constraint_rhs_storage,
            changing_constraint_rhs,
        )
    end

    added_variables = VariableRef[]
    added_constraints = ConstraintRef[]

    # Reformulate the objective function.
    _reformulate_objective(
        node, realizations,
        random_variables,
        added_variables,
        objective_storage,
        changing_objective_constant,
        changing_objective_coefficient,
    )

    # Reformulate fixed variables.
    for x in changing_variable_fixed_bound
        _reformulate_fixed_bound(
            node,
            realizations,
            random_variables,
            added_variables,
            added_constraints,
            variable_bound_storage,
            x,
        )
    end

    # Reformulate lower bounded variables.
    for x in changing_variable_lower_bound
        _reformulate_lower_bound(
            node,
            realizations,
            random_variables,
            added_variables,
            added_constraints,
            variable_bound_storage,
            x,
        )
    end

    # Reformulate upper bounded variables.
    for x in changing_variable_upper_bound
        _reformulate_upper_bound(
            node,
            realizations,
            random_variables,
            added_variables,
            added_constraints,
            variable_bound_storage,
            x,
        )
    end

    # Reformulate changing RHS term.
    for ci in changing_constraint_rhs
        _reformulate_constraint_rhs(
            node,
            realizations,
            random_variables,
            added_variables,
            constraint_rhs_storage,
            ci,
        )
    end

    return () -> begin
        delete(node.subproblem, added_variables)
        delete.(Ref(node.subproblem), added_constraints)
        set_objective_function(node.subproblem, node.stage_objective)
        return
    end
end

function _collect_changing_variable_bounds(
    node,
    variable_bound_storage,
    changing_variable_lower_bound,
    changing_variable_upper_bound,
    changing_variable_fixed_bound,
)
    bound = Dict{VariableRef, Any}()
    for x in all_variables(node.subproblem)
        l, u, f = -Inf, Inf, 0.0
        if has_lower_bound(x)
            l = lower_bound(x)
        end
        if has_upper_bound(x)
            u = upper_bound(x)
        end
        if is_fixed(x)
            f = fix_value(x)
        end
        if length(variable_bound_storage) >= 1
            if variable_bound_storage[1][x].l != l
                push!(changing_variable_lower_bound, x)
            end
            if variable_bound_storage[1][x].u != u
                push!(changing_variable_upper_bound, x)
            end
            if variable_bound_storage[1][x].f != f
                push!(changing_variable_fixed_bound, x)
            end
        end
        bound[x] = (l = l, u = u, f = f)
    end
    push!(variable_bound_storage, bound)
    return
end

function _collect_changing_objective(
    node,
    objective_storage,
    changing_objective_constant,
    changing_objective_coefficient,
)
    push!(objective_storage, convert(AffExpr, node.stage_objective))
    if length(objective_storage) > 1
        obj = objective_storage[end]
        if obj.constant != objective_storage[1].constant
            changing_objective_constant = true
        end
        for k in _dict_diff_keys(objective_storage[1].terms, obj.terms)
            push!(changing_objective_coefficient, k)
        end
    end
    return changing_objective_constant
end

function _collect_changing_constraint_rhs(
    node,
    constraint_rhs_storage,
    changing_constraint_rhs,
)
    rhs = Dict{ConstraintRef, Float64}()
    sets = (
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
        MOI.EqualTo{Float64}
    )
    for (F, S) in list_of_constraint_types(node.subproblem)
        if F == VariableRef || !(S in sets)
            continue
        end
        for ci in all_constraints(node.subproblem, F, S)
            obj = constraint_object(ci)
            rhs[ci] = MOI.constant(obj.set)
            if length(constraint_rhs_storage) >= 1
                if constraint_rhs_storage[1][ci] != rhs[ci]
                    push!(changing_constraint_rhs, ci)
                end
            end
        end
    end
    push!(constraint_rhs_storage, rhs)
    return
end

function _reformulate_objective(
    node::Node,
    realizations::Vector,
    random_variables::Vector{String},
    added_variables::Vector{VariableRef},
    objective_storage::Vector,
    changing_objective_constant::Bool,
    changing_objective_coefficient::Set{VariableRef},
)
    objective = convert(QuadExpr, copy(node.stage_objective))
    # Reformulate a changing objective constant.
    if changing_objective_constant
        new_name = "_SDDPjl_random_objective_constant_"
        y = _add_new_random_variable(
            node, new_name, random_variables, added_variables
        )
        for (r, o) in zip(realizations, objective_storage)
            r["support"][new_name] = o.constant
        end
        objective.aff.constant = 0.0
        objective.aff.terms[y] = 1.0
    end

    # No changes, so return current affine objective
    if length(changing_objective_coefficient) > 0
        # Reformulate changing objective coefficients.
        for x in changing_objective_coefficient
            new_name = "_SDDPjl_random_objective_$(name(x))_"
            y = _add_new_random_variable(
                node, new_name, random_variables, added_variables
            )
            for (r, o) in zip(realizations, objective_storage)
                r["support"][new_name] = get(o.terms, x, 0.0)
            end
            delete!.(Ref(objective.aff.terms), x)
            add_to_expression!(objective, 1.0, y, x)
        end
    end
    # Set the objective function to be written out.
    if length(objective.terms) == 0
        set_objective_function(node.subproblem, objective.aff)
    else
        set_objective_function(node.subproblem, objective)
    end
    return
end

function _reformulate_fixed_bound(
    ::Node,
    realizations::Vector,
    random_variables::Vector{String},
    ::Vector{VariableRef},
    ::Vector,
    variable_bound_storage::Vector,
    x::VariableRef,
)
    for (realization, bound) in zip(realizations, variable_bound_storage)
        realization["support"][name(x)] = bound[x].f
    end
    push!(random_variables, name(x))
    unfix(x)
end

function _reformulate_lower_bound(
    node::Node,
    realizations::Vector,
    random_variables::Vector{String},
    added_variables::Vector{VariableRef},
    added_constraints::Vector,
    variable_bound_storage::Vector,
    x::VariableRef,
)
    new_name = "_SDDPjl_lower_bound_$(name(x))_"
    y = _add_new_random_variable(
        node, new_name, random_variables, added_variables
    )
    c = @constraint(node.subproblem, x >= y)
    push!(added_constraints, c)
    delete_lower_bound(x)
    for (realization, bound) in zip(realizations, variable_bound_storage)
        realization["support"][new_name] = bound[x].l
    end
end

function _reformulate_upper_bound(
    node::Node,
    realizations::Vector,
    random_variables::Vector{String},
    added_variables::Vector{VariableRef},
    added_constraints::Vector,
    variable_bound_storage::Vector,
    x::VariableRef,
)
    new_name = "_SDDPjl_upper_bound_$(name(x))_"
    y = _add_new_random_variable(
        node, new_name, random_variables, added_variables
    )
    c = @constraint(node.subproblem, x <= y)
    push!(added_constraints, c)
    delete_upper_bound(x)
    for (realization, bound) in zip(realizations, variable_bound_storage)
        realization["support"][new_name] = bound[x].u
    end
end

function _reformulate_constraint_rhs(
    node,
    realizations,
    random_variables,
    added_variables,
    constraint_rhs_storage,
    ci,
)
    new_name = "_SDDPjl_rhs_$(name(ci))_"
    y = _add_new_random_variable(
        node, new_name, random_variables, added_variables
    )
    set_normalized_coefficient(ci, y, -1.0)
    set_normalized_rhs(ci, 0.0)
    for (realization, rhs) in zip(realizations, constraint_rhs_storage)
        realization["support"][new_name] = rhs[ci]
    end
    return
end

function _add_new_random_variable(
    node, new_name, random_variables, added_variables
)
    y = @variable(node.subproblem, base_name = new_name)
    push!(added_variables, y)
    push!(random_variables, new_name)
    return y
end


function _dict_diff_keys(
    x::AbstractDict{K, V}, y::AbstractDict{K, V}
) where {K, V}
    diff = Set{K}()
    for (k, v) in x
        if haskey(y, k)
            if v != y[k]
                push!(diff, k)
            end
        else
            push!(diff, k)
        end
    end
    for k in keys(y)
        if !haskey(x, k)
            push!(diff, k)
        end
    end
    return diff
end

function _subproblem_to_dict(subproblem::JuMP.Model)
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
    Base.read(
        io::IO,
        ::Type{PolicyGraph};
        bound::Float64 = 1e6,
    )::Tuple{PolicyGraph, TestScenarios}

Return a tuple containing a [`PolicyGraph`](@ref) object and a
[`TestScenarios`](@ref) read from `io` in the StochOptFormat file format.

See also: [`evaluate`](@ref).

WARNING: THIS FUNCTION IS EXPERIMENTAL. THINGS MAY CHANGE BETWEEN COMMITS. YOU
SHOULD NOT RELY ON THIS FUNCTIONALITY AS A LONG-TERM FILE FORMAT (YET).

In addition to potential changes to the underlying format, only a subset of
possible modifications are supported. These include:
- Additive random variables in the constraints or in the objective
- Multiplicative random variables in the objective

If your model uses something other than this, this function may throw an error
or silently build a non-convex model.

## Example

    open("my_model.sof.json", "r") do io
        model, test_scenarios = read(io, PolicyGraph)
    end
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
        if objective_sense(sp) != model_sense
            @warn(
                "Flipping the objective sense of node $(node_name) so that " *
                "it matches the majority of the subproblems."
            )
        end
        obj_sgn = objective_sense(sp) == model_sense ? 1 : -1
        objective_coefficients, objf = _convert_objective_function(
            sp,
            convert(
                Vector{String},
                data["nodes"][node_name]["random_variables"]
            )
        )
        parameterize(sp, Ω, P) do ω
            if ω !== nothing
                for (k, v) in ω
                    x = get(objective_coefficients, k, nothing)
                    if x !== nothing
                        if objf isa AffExpr
                            objf.terms[x.var] = x.aff + v * x.coef
                        else
                            objf.aff.terms[x.var] = x.aff + v * x.coef
                        end
                    end
                    fix(variable_by_name(sp, k), v)
                end
            end
            @stageobjective(sp, obj_sgn * objf)
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
    seekstart(io)
    SHA256 = bytes2hex(SHA.sha2_256(io))
    return model, _test_scenarios(data, SHA256)
end

function _test_scenarios(data::Dict, SHA256::String)
    substitute_nothing(x) = isempty(x) ? nothing : x
    scenarios = [
        TestScenario(
            scenario["probability"],
            [
                (item["node"], substitute_nothing(item["support"]))
                for item in scenario["scenario"]
            ]
        )
        for scenario in data["test_scenarios"]
    ]
    return TestScenarios(scenarios; SHA256 = SHA256)
end

function _convert_objective_function(sp::Model, rvs::Vector{String})
    return _convert_objective_function(sp, rvs, objective_function(sp))
end

function _convert_objective_function(sp::Model, ::Vector{String}, objf)
    return Dict{String, Any}(), objf
end

function _convert_objective_function(
    sp::Model, rvs::Vector{String}, objf::QuadExpr
)
    terms = Dict{String, Any}()
    aff_obj = copy(objf.aff)
    quad_terms = empty(copy(objf.terms))
    for (k, v) in objf.terms
        a, b = name(k.a), name(k.b)
        if a in rvs
            terms[a] = (var = k.b, coef = v, aff = get(aff_obj.terms, a, 0.0))
        elseif b in rvs
            terms[b] = (var = k.a, coef = v, aff = get(aff_obj.terms, b, 0.0))
        else
            quad_terms[k] = v
        end
    end
    if length(terms) == length(objf.terms)
        return terms, aff_obj
    else
        return terms, QuadExpr(aff_obj, quad_terms)
    end
end

"""
    write_to_file(
        model::PolicyGraph,
        filename::String;
        compression::MOI.FileFormats.AbstractCompressionScheme =
            MOI.FileFormats.AutomaticCompression(),
        kwargs...
    )

Write `model` to `filename` in the StochOptFormat file format.

Pass an argument to `compression` to override the default of automatically
detecting the file compression to use based on the extension of `filename`.

See [`Base.write(::IO, ::PolicyGraph)`](@ref) for information on the
keyword arguments that can be provided.

WARNING: THIS FUNCTION IS EXPERIMENTAL. SEE THE FULL WARNING IN
[`Base.write(::IO, ::PolicyGraph)`](@ref).

## Example

    write_to_file(model, "my_model.sof.json"; test_scenarios = 10)
"""
function write_to_file(
    model::PolicyGraph,
    filename::String;
    compression::MOI.FileFormats.AbstractCompressionScheme =
        MOI.FileFormats.AutomaticCompression(),
    kwargs...
)
    return MOI.FileFormats.compressed_open(filename, "w", compression) do io
        Base.write(io, model; kwargs...)
    end
end

"""
    read_from_file(
        filename::String;
        compression::MOI.FileFormats.AbstractCompressionScheme =
            MOI.FileFormats.AutomaticCompression(),
        kwargs...
    )::Tuple{PolicyGraph, TestScenarios}

Return a tuple containing a [`PolicyGraph`](@ref) object and a
[`TestScenarios`](@ref) read from `filename` in the StochOptFormat file format.

Pass an argument to `compression` to override the default of automatically
detecting the file compression to use based on the extension of `filename`.

See [`Base.read(::IO, ::Type{PolicyGraph})`](@ref) for information on the
keyword arguments that can be provided.

WARNING: THIS FUNCTION IS EXPERIMENTAL. SEE THE FULL WARNING IN
[`Base.read(::IO, ::Type{PolicyGraph})`](@ref).

## Example

    model, test_scenarios = read_from_file("my_model.sof.json")
"""
function read_from_file(
    filename::String;
    compression::MOI.FileFormats.AbstractCompressionScheme =
        MOI.FileFormats.AutomaticCompression(),
    kwargs...
)
    return MOI.FileFormats.compressed_open(filename, "r", compression) do io
        Base.read(io, PolicyGraph; kwargs...)
    end
end

"""
    evaluate(
        model::PolicyGraph{T}, test_scenarios::TestScenarios{T, S}
    ) where {T, S}

Evaluate the performance of the policy contained in `model` after a call to
[`train`](@ref) on the scenarios specified by `test_scenarios`.

## Example

    model, test_scenarios = read_from_file("my_model.sof.json")
    train(model; iteration_limit = 100)
    simulations = evaluate(model, test_scenarios)
"""
function evaluate(
    model::PolicyGraph{T}, test_scenarios::TestScenarios{T, S}
) where {T, S}
    test_scenarios.last = 0
    simulations = simulate(
        model,
        length(test_scenarios.scenarios);
        sampling_scheme = test_scenarios,
        custom_recorders = Dict{Symbol, Function}(
            :primal => (sp) -> begin
                Dict{String, Float64}(
                    name(x) => value(x) for x in all_variables(sp)
                    if !isempty(name(x))
                )
            end
        )
    )
    return Dict(
        "problem_sha256_checksum" => test_scenarios.SHA256,
        "scenarios" => [
            [
                Dict{String, Any}(
                    "objective" => s[:stage_objective],
                    "primal" => s[:primal]
                )
                for s in sim
            ] for sim in simulations
        ]
    )
end
