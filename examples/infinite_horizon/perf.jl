using Kokako, Gurobi

function perf()
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
            bellman_function = Kokako.AverageCut(lower_bound = 0.0),
            optimizer = with_optimizer(Gurobi.Optimizer, OutputFlag = 0),
            sense = :Min
            ) do subproblem, node
        @variable(subproblem, 0 <= r <= 2, Kokako.State, (initial_value = 1))
        @constraint(subproblem, r.out - r.in == 0.1)
        @stageobjective(subproblem, 100 * r.out)
    end
    Kokako.train(model, iteration_limit = 1)
    return model
end


function get_outgoing_state_1(node::Kokako.Node)
    values = Dict{Symbol, Float64}()
    for (name, state) in node.states
        outgoing_value = JuMP.value(state.out)
        if JuMP.has_upper_bound(state.out)
            current_bound = JuMP.upper_bound(state.out)
            if current_bound < outgoing_value
                outgoing_value = current_bound
            end
        elseif JuMP.has_lower_bound(state.out)
            current_bound = JuMP.lower_bound(state.out)
            if current_bound > outgoing_value
                outgoing_value = current_bound
            end
        end
        values[name] = outgoing_value
    end
    return values
end

function get_outgoing_state_2(node::Kokako.Node)
    values = Dict{Symbol, Float64}()
    for (name, state) in node.states
        values[name] = JuMP.value(state.out)
    end
    return values
end

model = perf()

using BenchmarkTools

# @btime get_outgoing_state_1($(model.nodes[1]))
# @btime get_outgoing_state_2($(model.nodes[1]))

x = model.nodes[1].states[:r].out

# @btime JuMP.upper_bound($x)
# @btime JuMP.has_upper_bound($x)
# @btime JuMP.lower_bound($x)
# @btime JuMP.has_lower_bound($x)
@btime JuMP.value($x)

foo(m, x) = MOI.get(m, MOI.VariablePrimal(), x)
@btime foo($(model.nodes[1].subproblem), $x)

# # Internal function: get the values of the dual variables associated with the
# # fixed incoming state variables. Requires node.subproblem to have been solved
# # with DualStatus == FeasiblePoint.
# function get_dual_variables(node::Node)
#     # Note: due to JuMP's dual convention, we need to flip the sign for
#     # maximization problems.
#     dual_sign = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1.0 : -1.0
#     values = Dict{Symbol, Float64}()
#     for (name, state) in node.states
#         ref = JuMP.FixRef(state.in)
#         values[name] = dual_sign * JuMP.dual(ref)
#     end
#     return values
# end
