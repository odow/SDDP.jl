#  Copyright 2017-20, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# ============================================================================ #
#
#   Code to implement a JuMP variable extension.
#
#   Usage:
#   julia> @variable(subproblem, 0 <= x[i=1:2] <= i,
#              SDDP.State, initial_value = i)
#
#   julia> x
#   2-element Array{State{VariableRef},1}:
#     State(x[1]_in,x[1]_out)
#     State(x[2]_in,x[2]_out)
#
#   julia> x[1].in
#   x[1]_in
#
#   julia> typeof(x[1].in)
#   VariableRef
#
#   julia> x[2].out
#   x[2]_out
#
#   Assuming subproblem has been solved, and there exists a primal solution
#   julia> x_values = JuMP.value.(x)
#   2-element Array{State{Float64},1}:
#     State(0.0,1.0)
#     State(1.2,3.0)
#
#   julia> x_values[1].out
#   1.0
# ============================================================================ #

struct StateInfo
    in::JuMP.VariableInfo
    out::JuMP.VariableInfo
    initial_value::Float64
    kwargs
end

function JuMP.build_variable(
    _error::Function,
    info::JuMP.VariableInfo,
    ::Type{State};
    initial_value = NaN,
    kwargs...,
)
    if isnan(initial_value)
        _error(
            "When creating a state variable, you must set the " *
            "`initial_value` keyword to the value of the state variable at" *
            " the root node.",
        )
    end
    return StateInfo(
        JuMP.VariableInfo(
            false,
            NaN,  # lower bound
            false,
            NaN,  # upper bound
            false,
            NaN,  # fixed value
            false,
            NaN,  # start value
            false,
            false, # binary and integer
        ),
        info,
        initial_value,
        kwargs,
    )
end

function JuMP.add_variable(subproblem::JuMP.Model, state_info::StateInfo, name::String)
    state = State(
        JuMP.add_variable(subproblem, JuMP.ScalarVariable(state_info.in), name * "_in"),
        JuMP.add_variable(subproblem, JuMP.ScalarVariable(state_info.out), name * "_out"),
    )
    integrality_handler = get_integrality_handler(subproblem)
    setup_state(subproblem, state, state_info, name, integrality_handler)
    return state
end

JuMP.variable_type(model::JuMP.Model, ::Type{State}) = State

function JuMP.value(state::State{JuMP.VariableRef})
    return State(JuMP.value(state.in), JuMP.value(state.out))
end

# Overload for broadcast syntax such as `JuMP.value.([state_1, state_2])`.
Broadcast.broadcastable(state::State{JuMP.VariableRef}) = Ref(state)

# ==============================================================================

function JuMP.set_optimizer(model::SDDP.PolicyGraph, optimizer)
    for node in values(model.nodes)
        set_optimizer(node.subproblem, optimizer)
    end
end
