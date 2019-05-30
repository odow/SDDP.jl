#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

struct ScenarioTreeNode{T}
    node::Node{T}
    noise::Any
    probability::Float64
    children::Vector{ScenarioTreeNode{T}}
    states::Dict{Symbol, State{JuMP.VariableRef}}
end

struct ScenarioTree{T}
    children::Vector{ScenarioTreeNode{T}}
end

function add_node_to_scenario_tree(
    parent::Vector{ScenarioTreeNode{T}}, pg::PolicyGraph{T},
    node::Node{T}, probability::Float64, added_nodes::Set{T}
) where {T}
    if node.objective_state !== nothing
        error("Objective states detected! Unable to formulate deterministic equivalent.")
    elseif node.belief_state !== nothing
        error("Belief states detected! Unable to formulate deterministic equivalent.")
    end
    for noise in node.noise_terms
        scenario_node = ScenarioTreeNode(
            node,
            noise.term,
            probability * noise.probability,
            ScenarioTreeNode{T}[],
            Dict{Symbol, State{JuMP.VariableRef}}()
        )
        for child in node.children
            if child.term in added_nodes
                error(
                    "Cycle detected in the policy graph! Unable to formulate " *
                    "deterministic equivalent."
                )
            end
            push!(added_nodes, child.term)
            add_node_to_scenario_tree(
                scenario_node.children, pg, pg[child.term],
                probability * noise.probability * child.probability, added_nodes
            )
        end
        push!(parent, scenario_node)
    end
    return
end

function copy_and_replace_variables(
    src::JuMP.VariableRef,
    src_to_dest_variable::Dict{JuMP.VariableRef, JuMP.VariableRef}
)
    return src_to_dest_variable[src]
end

function copy_and_replace_variables(
    src::JuMP.GenericAffExpr,
    src_to_dest_variable::Dict{JuMP.VariableRef, JuMP.VariableRef}
)
    return JuMP.GenericAffExpr(
        src.constant,
        (src_to_dest_variable[key] => val for (key, val) in src.terms)...
    )
end

function add_scenario_to_ef(model::JuMP.Model, child::ScenarioTreeNode)
    node = child.node
    parameterize(node, child.noise)
    # Add variables:
    src_variables = JuMP.all_variables(node.subproblem)
    x = @variable(model, [1:length(src_variables)])
    # TODO: set nice names.
    var_src_to_dest = Dict(s => d for (s, d) in zip(src_variables, x))
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
    subproblem_objective = copy_and_replace_variables(
        node.stage_objective, var_src_to_dest
    )
    JuMP.set_objective_function(
        model, current + child.probability * subproblem_objective
    )
    # Add state variables to child.states:
    for (key, state) in node.states
        child.states[key] = State(
            var_src_to_dest[state.in],
            var_src_to_dest[state.out]
        )
    end
    # Recurse down the tree.
    for child_2 in child.children
        add_scenario_to_ef(model, child_2)
    end
    return
end

function add_linking_constraints(model::JuMP.Model, node::ScenarioTreeNode)
    for child in node.children
        for key in keys(node.states)
            @constraint(model, node.states[key].out == child.states[key].in)
        end
        add_linking_constraints(model, child)
    end
end

function deterministic_equivalent(pg::PolicyGraph{T}, optimizer=nothing) where {T}
    # Step 1: convert the policy graph into a scenario tree.
    tree = ScenarioTree{T}(ScenarioTreeNode{T}[])
    added_nodes = Set{T}()
    for child in pg.root_children
        add_node_to_scenario_tree(
            tree.children, pg, pg[child.term], child.probability, added_nodes
        )
    end
    # Step 2: create a extensive-form JuMP model and add subproblems.
    model = optimizer === nothing ? JuMP.Model() : JuMP.Model(optimizer)
    for child in tree.children
        add_scenario_to_ef(model, child)
    end
    # Step 3: add linking constraints between the nodes in the scenario tree.
    for child in tree.children
        add_linking_constraints(model, child)
        for (key, value) in pg.initial_root_state
            JuMP.fix(child.states[key].in, value; force=true)
        end
    end
    # Return the model
    return model
end
