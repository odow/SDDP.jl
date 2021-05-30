#  Copyright 2017-21, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function throw_detequiv_error(msg::String)
    return error("Unable to formulate deterministic equivalent: ", msg)
end

struct ScenarioTreeNode{T}
    node::Node{T}
    noise::Any
    probability::Float64
    children::Vector{ScenarioTreeNode{T}}
    states::Dict{Symbol,State{JuMP.VariableRef}}
end

struct ScenarioTree{T}
    children::Vector{ScenarioTreeNode{T}}
end

function add_node_to_scenario_tree(
    parent::Vector{ScenarioTreeNode{T}},
    pg::PolicyGraph{T},
    node::Node{T},
    probability::Float64,
    check_time_limit::Function,
) where {T}
    if node.objective_state !== nothing
        throw_detequiv_error("Objective states detected!")
    elseif node.belief_state !== nothing
        throw_detequiv_error("Belief states detected!")
    elseif length(node.bellman_function.global_theta.cut_oracle.cuts) > 0
        throw_detequiv_error(
            "Model has been used for training. Can only form deterministic " *
            "equivalent on a fresh model.",
        )
    else
        check_time_limit()
    end
    for noise in node.noise_terms
        scenario_node = ScenarioTreeNode(
            node,
            noise.term,
            probability * noise.probability,
            ScenarioTreeNode{T}[],
            Dict{Symbol,State{JuMP.VariableRef}}(),
        )
        for child in node.children
            add_node_to_scenario_tree(
                scenario_node.children,
                pg,
                pg[child.term],
                probability * noise.probability * child.probability,
                check_time_limit,
            )
        end
        push!(parent, scenario_node)
    end
    return
end

function copy_and_replace_variables(
    src::Vector,
    map::Dict{JuMP.VariableRef,JuMP.VariableRef},
)
    return copy_and_replace_variables.(src, Ref(map))
end

function copy_and_replace_variables(
    src::Real,
    ::Dict{JuMP.VariableRef,JuMP.VariableRef},
)
    return src
end

function copy_and_replace_variables(
    src::JuMP.VariableRef,
    src_to_dest_variable::Dict{JuMP.VariableRef,JuMP.VariableRef},
)
    return src_to_dest_variable[src]
end

function copy_and_replace_variables(
    src::JuMP.GenericAffExpr,
    src_to_dest_variable::Dict{JuMP.VariableRef,JuMP.VariableRef},
)
    return JuMP.GenericAffExpr(
        src.constant,
        Pair{VariableRef,Float64}[
            src_to_dest_variable[key] => val for (key, val) in src.terms
        ],
    )
end

function copy_and_replace_variables(
    src::JuMP.GenericQuadExpr,
    src_to_dest_variable::Dict{JuMP.VariableRef,JuMP.VariableRef},
)
    return JuMP.GenericQuadExpr(
        copy_and_replace_variables(src.aff, src_to_dest_variable),
        Pair{UnorderedPair{VariableRef},Float64}[
            UnorderedPair{VariableRef}(
                src_to_dest_variable[pair.a],
                src_to_dest_variable[pair.b],
            ) => coef for (pair, coef) in src.terms
        ],
    )
end

function copy_and_replace_variables(
    src::Any,
    ::Dict{JuMP.VariableRef,JuMP.VariableRef},
)
    return throw_detequiv_error(
        "`copy_and_replace_variables` is not implemented for functions like `$(src)`.",
    )
end

function add_scenario_to_ef(
    model::JuMP.Model,
    child::ScenarioTreeNode,
    check_time_limit::Function,
)
    check_time_limit()
    node = child.node
    parameterize(node, child.noise)
    # Add variables:
    src_variables = JuMP.all_variables(node.subproblem)
    x = @variable(model, [1:length(src_variables)])
    var_src_to_dest = Dict{JuMP.VariableRef,JuMP.VariableRef}()
    for (src, dest) in zip(src_variables, x)
        var_src_to_dest[src] = dest
        JuMP.set_name(dest, JuMP.name(src))
    end
    # Add constraints:
    for (F, S) in JuMP.list_of_constraint_types(node.subproblem)
        for con in JuMP.all_constraints(node.subproblem, F, S)
            obj = JuMP.constraint_object(con)
            new_func = copy_and_replace_variables(obj.func, var_src_to_dest)
            @constraint(model, new_func in obj.set)
        end
    end
    # Add objective:
    current = JuMP.objective_function(model)
    subproblem_objective =
        copy_and_replace_variables(node.stage_objective, var_src_to_dest)
    JuMP.set_objective_function(
        model,
        current + child.probability * subproblem_objective,
    )
    # Add state variables to child.states:
    for (key, state) in node.states
        child.states[key] =
            State(var_src_to_dest[state.in], var_src_to_dest[state.out])
    end
    # Recurse down the tree.
    for child_2 in child.children
        add_scenario_to_ef(model, child_2, check_time_limit)
    end
    return
end

function add_linking_constraints(
    model::JuMP.Model,
    node::ScenarioTreeNode,
    check_time_limit::Function,
)
    check_time_limit()
    for child in node.children
        for key in keys(node.states)
            @constraint(model, node.states[key].out == child.states[key].in)
        end
        add_linking_constraints(model, child, check_time_limit)
    end
end

"""
    deterministic_equivalent(
        pg::PolicyGraph{T},
        optimizer = nothing;
        time_limit::Union{Real, Nothing} = 60.0
    )

Form a JuMP model that represents the deterministic equivalent of the problem.

## Examples

    deterministic_equivalent(model)
    deterministic_equivalent(model, GLPK.Optimizer)
"""
function deterministic_equivalent(
    pg::PolicyGraph{T},
    optimizer = nothing;
    time_limit::Union{Real,Nothing} = 60.0,
) where {T}
    # Step 0: helper function for the time limit.
    start_time = time()
    time_limit = time_limit === nothing ? typemax(Float64) : Float64(time_limit)
    function check_time_limit()
        if time() - start_time > time_limit::Float64
            throw_detequiv_error("Time limit exceeded!")
        end
    end
    # Step 1: convert the policy graph into a scenario tree.
    if is_cyclic(pg)
        throw_detequiv_error("Cyclic policy graph detected!")
    end
    tree = ScenarioTree{T}(ScenarioTreeNode{T}[])
    for child in pg.root_children
        add_node_to_scenario_tree(
            tree.children,
            pg,
            pg[child.term],
            child.probability,
            check_time_limit,
        )
    end
    # Step 2: create a extensive-form JuMP model and add subproblems.
    model = optimizer === nothing ? JuMP.Model() : JuMP.Model(optimizer)
    for child in tree.children
        add_scenario_to_ef(model, child, check_time_limit)
    end
    # Step 3: add linking constraints between the nodes in the scenario tree.
    for child in tree.children
        add_linking_constraints(model, child, check_time_limit)
        for (key, value) in pg.initial_root_state
            JuMP.fix(child.states[key].in, value; force = true)
        end
    end
    # Return the model
    return model
end
