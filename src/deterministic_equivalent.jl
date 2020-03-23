#  Copyright 2017-20, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function throw_detequiv_error(msg::String)
    error("Unable to formulate deterministic equivalent: ", msg)
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

function is_cyclic(G::PolicyGraph{T}) where {T}
    # We implement Tarjan's strongly connected components algorithm to detect
    # cycles in a directed graph in O(|V| + |E|) time. See this Wiki for details
    # https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
    # The notation here follows the pseudocode in the Wikipedia article, rather
    # than the typical JuMP style guide.
    #
    # Since we're only checking for cyclic graphs, we can stop as soon as on is
    # found. A cyclic graph has a stongly connected component with at least two
    # components, or it has a node with connects to itself. That means we don't
    # need to store the set of all strongly connected components.
    index_counter = 0
    S = T[]
    low_link = Dict{T,Int}()
    index = Dict{T,Int}()
    on_stack = Dict{T,Bool}()
    function strong_connect(v)
        index[v] = index_counter
        low_link[v] = index_counter
        index_counter += 1
        push!(S, v)
        on_stack[v] = true
        for child in G[v].children
            w = child.term
            if v == w
                # Cycle detected: Type I: a node that loops to itself.
                return true
            end
            if !haskey(index, w)
                if strong_connect(w)
                    # A cycle was detected further down the tree. Propogate it
                    # upwards.
                    return true
                end
                low_link[v] = min(low_link[v], low_link[w])
            elseif on_stack[w]
                low_link[v] = min(low_link[v], index[w])
            end
        end
        if low_link[v] == index[v]
            scc = T[]
            w = G.root_node
            while v != w
                w = pop!(S)
                on_stack[w] = false
                push!(scc, w)
            end
            if length(scc) > 1
                # Cycle detected: Type II: a strongly connected component with
                # more than one element.
                return true
            end
        end
        return false  # No cycle detected.
    end
    for v in keys(G.nodes)
        if !haskey(index, v)
            if strong_connect(v)
                # Cycle detected!
                return true
            end
        end
    end
    return false
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
            "equivalent on a fresh model."
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

copy_and_replace_variables(src::Real, ::Dict{JuMP.VariableRef,JuMP.VariableRef}) = src

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

function copy_and_replace_variables(src::Any, ::Dict{JuMP.VariableRef,JuMP.VariableRef})
    throw_detequiv_error("`copy_and_replace_variables` is not implemented for functions like `$(src)`.")
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
    subproblem_objective = copy_and_replace_variables(node.stage_objective, var_src_to_dest)
    JuMP.set_objective_function(model, current + child.probability * subproblem_objective)
    # Add state variables to child.states:
    for (key, state) in node.states
        child.states[key] = State(var_src_to_dest[state.in], var_src_to_dest[state.out])
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
