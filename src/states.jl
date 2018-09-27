function add_state_variable(subproblem::JuMP.Model,
                            name::Symbol,
                            incoming::JuMP.VariableRef,
                            outgoing::JuMP.VariableRef)
    states = extension(subproblem).states
    JuMP.fix(incoming, 0.0)  # Fix the variable value.
    if haskey(states, name)
        error("The state $(name) already exists.")
    else
        states[name] = State(incoming, outgoing)
    end
    return
end
