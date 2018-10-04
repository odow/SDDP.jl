struct StateInfo
    in::JuMP.VariableInfo
    out::JuMP.VariableInfo
    root_value::Float64
end

function JuMP.build_variable(
    _error::Function, info::JuMP.VariableInfo, ::Type{State};
    root_value = NaN, kwargs...)
    if isnan(root_value)
        _error("When creating a state variable, you must set the `root_value`" *
               " keyword to the value of the state variable at the root node.")
    end
    return StateInfo(
        JuMP.VariableInfo(
            false, NaN,  # lower bound
            false, NaN,  # upper bound
            false, NaN,  # fixed value
            false, NaN,  # start value
            false, false # binary and integer
        ),
        info,
        root_value
    )
end

function JuMP.add_variable(
        subproblem::JuMP.Model, state_info::StateInfo, name::String)
    state = State(
        JuMP.add_variable(
            subproblem, JuMP.ScalarVariable(state_info.in), name * "_in"),
        JuMP.add_variable(
            subproblem, JuMP.ScalarVariable(state_info.out), name * "_out")
    )
    node = get_node(subproblem)
    sym_name = Symbol(name)
    if haskey(node.states, sym_name)
        error("The state $(sym_name) already exists.")
    end
    node.states[sym_name] = state
    graph = get_policy_graph(subproblem)
    graph.initial_root_state[sym_name] = state_info.root_value
    return state
end

JuMP.variable_type(model::JuMP.Model, ::Type{State}) = State

function JuMP.result_value(state::State{JuMP.VariableRef})
    return State(JuMP.result_value(state.in), JuMP.result_value(state.out))
end

# using JuMP, GLPK
# model = Model(with_optimizer(GLPK.Optimizer))
# @variable(model, 0 <= x[i=1:3] <= i, State, root_value = i)
# JuMP.optimize!(model)
# JuMP.result_value.(x)
