#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Revise

using SDDP
import Distributions
import HiGHS
import Gurobi
import Random
import Plots

function add_ddu_matrices(sp, matrices::Vector{<:Matrix}; M)
    N = length(matrices)
    @variable(sp, __ddu__[1:N], Bin)
    @constraint(sp, sum(__ddu__) == 1)
    sp.ext[:__ddu__] = (M = M, y = __ddu__, matrices = matrices)
    return __ddu__
end

function solve_decision_dependent_trajectory(
    model::SDDP.PolicyGraph{Int},
    incoming_state_value,
    variables::Vector{Symbol} = Symbol[];
    explore::Bool = true,
    depth::Union{Nothing,Int} = nothing,
)
    for i in keys(model.nodes)
        model[i].ext[:_ddu_is_set] = true
    end
    function sample_node(Φ::Matrix{<:Real}, y::Int)
        r = rand()
        for j in 1:size(Φ, 2)
            r -= Φ[y, j]
            if r <= 0
                return j
            end
        end
        return nothing
    end
    sampled_states = Dict{Symbol,Float64}[]
    belief_states = Tuple{Int,Dict{Int,Float64}}[]
    current_belief = SDDP.initialize_belief(model)
    cumulative_value = 0.0
    scenario_path = Tuple{Int,Any}[]
    simulation = Dict{Symbol,Any}[]
    i, y = 0, 1
    Φ = first(values(model.nodes)).subproblem.ext[:__ddu__].matrices
    while length(scenario_path) < something(depth, Inf)
        if depth === nothing
            i = sample_node(Φ[y], i + 1)
            if i === nothing
                break
            end
        else
            j = nothing
            while j === nothing
                j = sample_node(Φ[y], i + 1)
            end
            i = j
        end
        node = model[i]
        ω = SDDP.sample_noise(node.noise_terms)
        push!(scenario_path, (i, ω))
        if node.belief_state !== nothing
            belief = node.belief_state::SDDP.BeliefState{Int}
            current_belief = belief.updater(
                node,
                belief.belief,
                current_belief,
                belief.partition_index,
                ω,
            )
            push!(belief_states, (belief.partition_index, copy(current_belief)))
        end
        subproblem_results = SDDP.solve_subproblem(
            model,
            node,
            incoming_state_value,
            ω,
            scenario_path,
            duality_handler = nothing,
        )
        __ddu__ = node.subproblem.ext[:__ddu__]
        y = findfirst([round(Bool, value(y)) for y in __ddu__.y])
        if explore && rand() < 0.8
            y = rand(1:length(__ddu__.y))
        end
        cumulative_value += subproblem_results.stage_objective
        incoming_state_value = copy(subproblem_results.state)
        push!(sampled_states, incoming_state_value)
        # Record useful variables from the solve.
        store = Dict{Symbol,Any}(
            :node_index => i,
            :noise_term => ω,
            :stage_objective => subproblem_results.stage_objective,
            :bellman_term =>
                subproblem_results.objective -
                subproblem_results.stage_objective,
            # :objective_state => objective_state_vector,
            :belief => copy(current_belief),
        )
        # Loop through the primal variable values that the user wants.
        for variable in variables
            if haskey(node.subproblem.obj_dict, variable)
                # Note: we broadcast the call to value for variables which are
                # containers (like Array, Containers.DenseAxisArray, etc). If
                # the variable is a scalar (e.g. just a plain VariableRef), the
                # broadcast preseves the scalar shape.
                # TODO: what if the variable container is a dictionary? They
                # should be using Containers.SparseAxisArray, but this might not
                # always be the case...
                store[variable] = JuMP.value.(node.subproblem[variable])
            elseif skip_undefined_variables
                store[variable] = NaN
            else
                error(
                    "No variable named $(variable) exists in the subproblem.",
                    " If you want to simulate the value of a variable, make ",
                    "sure it is defined in _all_ subproblems, or pass ",
                    "`skip_undefined_variables=true` to `simulate`.",
                )
            end
        end
        push!(simulation, store)
        if depth === nothing &&
           rand() >= sum(child.probability for child in node.children; init = 0)
            break
        end
    end
    return (
        scenario_path = scenario_path,
        sampled_states = sampled_states,
        objective_states = NTuple{0,Float64}[],
        belief_states = belief_states,
        cumulative_value = cumulative_value,
        simulation = simulation,
    )
end

struct DecisionDependentForwardPass <: SDDP.AbstractForwardPass end

function SDDP.forward_pass(
    model::SDDP.PolicyGraph{Int},
    options::SDDP.Options,
    ::DecisionDependentForwardPass,
)
    incoming_state_value = copy(options.initial_state)
    return solve_decision_dependent_trajectory(model, incoming_state_value)
end


Φ(ρ, ε, z) = [
    #= R   =# 0.5 0.5 0 0
    #= q_L =# ρ*(1-ε)*(1-z) ρ*ε*(1-z) z 0
    #= q_H =# ρ*ε*(1-z) ρ*(1-ε)*(1-z) 0 z
    #= d_L =# ρ*(1-ε) ρ*ε 0 0
    #= d_H =# ρ*ε ρ*(1-ε) 0 0
]

ρ, ε = 0.9, 0.0
c_market, c_price = 3, 1.0
Ω = [
    [0, 2, 4, 6, 8, 10],
    [0, 2, 4, 6, 8, 10],
    [5, 10],
    [5, 10],
]
P = [
    fill(1 / length(Ω[1]), length(Ω[1])),
    fill(1 / length(Ω[2]), length(Ω[2])),
    [0.8, 0.2],
    [0.2, 0.8]
]
graph = SDDP.Graph(0)
SDDP.add_node.((graph,), 1:4)
SDDP.add_ambiguity_set(graph, [1, 2], 1e3)
SDDP.add_ambiguity_set(graph, [3, 4], 1e3)
Φ̅ = Φ(ρ, ε, 0.5)
for i in 1:5, j in 1:4
    Φ̅[i, j] > 0 && SDDP.add_edge(graph, (i-1) => j, Φ̅[i, j])
end
model = SDDP.PolicyGraph(
    graph;
    sense = :Max,
    optimizer = Gurobi.Optimizer,
    upper_bound = c_price * maximum(Ω[4]) / (1 - ρ),
) do sp, node
    @variable(sp, x >= 0, SDDP.State, initial_value = 0)
    @variable(sp, u_sell >= 0)
    z = add_ddu_matrices(sp, [Φ(ρ, ε, 0), Φ(ρ, ε, 1)]; M = 1e4)
    sp[:z] = z
    @constraint(sp, con_balance, x.out == x.in - u_sell + 0.0)
    if node in (1, 2)  # farm
        @stageobjective(sp, -c_market * z[2])
        SDDP.parameterize(sp, Ω[node], P[node]) do ω
            return set_normalized_rhs(con_balance, ω)
        end
    else  # market
        @stageobjective(sp, c_price * u_sell)
        SDDP.parameterize(sp, Ω[node], P[node]) do ω
            return set_upper_bound(u_sell, ω)
        end
    end
end
SDDP.train(
    model;
    duality_handler = SDDP.LagrangianDuality(),
    cut_type = SDDP.MULTI_CUT,
    forward_pass = DecisionDependentForwardPass(),
    iteration_limit = 10,
    log_every_iteration = true,
)

ret = solve_decision_dependent_trajectory(
    model,
    model.initial_root_state,
    [:x, :u_sell, :z];
    explore = false,
    depth = 50,
)

Plots.plot(
    Plots.plot(
        map(d -> d[:node_index], ret.simulation);
        ylabel = "Node",
        linetype = :step,
    ),
    Plots.plot(
        map(d -> (d[:belief][1] + d[:belief][3]), ret.simulation);
        ylabel = "Belief(low)",
        linetype = :step,
        ylims = (0, 1),
    ),
    Plots.plot(
        map(d -> d[:x].out, ret.simulation);
        ylabel = "x",
        ylims = (0, maximum(d -> d[:x].out, ret.simulation)),
        linetype = :step,
    ),
    Plots.plot(
        map(d -> d[:z][2], ret.simulation),
        ylims = (0, 1),
        ylabel = "z",
        linetype = :step,
    ),
    Plots.plot(
        map(d -> d[:u_sell], ret.simulation);
        ylabel = "u_sell",
        linetype = :step,
    ),
    xlims = (0, length(ret.simulation) + 1),
    legend = false,
    xlabel = "Simulation step",
)

args = NTuple{3,Float64}[
    (data[:x].out, data[:belief][1] + data[:belief][3], data[:z][2])
    for data in ret.simulation if data[:node_index] in (1, 2)
]
Plots.scatter(args; xlabel = "x", ylabel = "z", legend = false)

function plot_simulation(ret)
    x = Tuple{Float64,Float64}[]
    u = Tuple{Float64,Float64}[]
    q = Tuple{Float64,Float64}[]
    d = Tuple{Float64,Float64}[]
    t_index = 1.0
    for data in ret.simulation
        if data[:node_index] in (1, 2)
            t_index += 1.0
            push!(x, (t_index, data[:x].out))
            push!(q, (t_index, data[:noise_term]))
        else
            push!(u, (t_index+0.1, -data[:u_sell]))
            push!(d, (t_index-0.1, -data[:noise_term]))
        end
    end
    plot_x = Plots.plot(x)
    plot_u = Plots.bar(
        u;
        xlabel = "Week",
        ylabel = "Quantity [kg]",
        ylims = (-maximum(Ω[3]) - 1, maximum(Ω[1]) + 1),
        bar_width = 0.2,
        label = "u_sell",
    )
    Plots.bar!(plot_u, q; bar_width = 0.4, label = "Supply")
    Plots.bar!(plot_u, d; bar_width = 0.2, label = "Demand")
    return Plots.plot(plot_x, plot_u)
end

# plot_simulation(ret)
