#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Revise

using SDDP
import Distributions
import HiGHS
import Random

function _add_ddu_constraints(model::SDDP.PolicyGraph{Int}, i::Int)
    node = model[i]
    if get(node.ext, :_ddu_is_set, false)
        return
    end
    nominal_P = [
        child.probability * noise.probability
        for child in node.children for noise in model[child.term].noise_terms
    ]
    push!(node.bellman_function.risk_set_cuts, nominal_P)
    N = length(nominal_P)
    SDDP._add_locals_if_necessary(node, node.bellman_function, N)
    θʲʷ = [node.bellman_function.local_thetas[i].theta for i in 1:N]
    Θ = node.bellman_function.global_theta.theta
    ddu = node.subproblem.ext[:__ddu__]
    for (d, y_d) in enumerate(ddu.y)
        P_d = [
            ddu.matrices[d][i+1, child.term] * noise.probability
            for child in node.children
            for noise in model[child.term].noise_terms
        ]
        slack = ddu.M * (1 - y_d)
        if JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE
            JuMP.@constraint(node.subproblem, Θ >= P_d' * θʲʷ - slack)
        else
            JuMP.@constraint(node.subproblem, Θ <= P_d' * θʲʷ + slack)
        end
    end
    node.ext[:_ddu_is_set] = true
    return
end

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
)
    for i in keys(model.nodes)
        _add_ddu_constraints(model, i)
    end
    function sample_node(Φ::Matrix{Float64}, y::Int)
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
    cumulative_value = 0.0
    scenario_path = Tuple{Int,Any}[]
    simulation = Dict{Symbol,Any}[]
    i, y = 0, 1
    Φ = first(values(model.nodes)).subproblem.ext[:__ddu__].matrices
    while true
        i = sample_node(Φ[y], i + 1)
        if i === nothing
            break
        end
        node = model[i]
        ω = SDDP.sample_noise(node.noise_terms)
        push!(scenario_path, (i, ω))
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
            # :belief => copy(current_belief),
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
        if rand() <= 1 - sum(child.probability for child in node.children)
            break
        end
    end
    return (
        scenario_path = scenario_path,
        sampled_states = sampled_states,
        objective_states = NTuple{0,Float64}[],
        belief_states = Tuple{Int,Dict{Int,Float64}}[],
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

function build_and_cheese_model(seed::Int)
    Φ(z, ρ) = [1 0; ρ * (1 - z) z; ρ 0]
    ρ = 0.96
    graph = SDDP.Graph(0)
    SDDP.add_node.((graph,), 1:2)
    Φ̅ = Φ(0.5, ρ)
    for i in 1:3, j in 1:2
        Φ̅[i, j] > 0 && SDDP.add_edge(graph, (i-1) => j, Φ̅[i, j])
    end
    Random.seed!(seed)
    # Ω_q = sort!(rand(Distributions.Uniform(0, 20), 2))
    # Ω_d = sort!(rand(Distributions.Uniform(10, 50), 3))
    return SDDP.PolicyGraph(
        graph;
        sense = :Max,
        upper_bound = 1e3, # 50 * 2 / (1 - ρ),
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, 0 <= x <= 20, SDDP.State, initial_value = 0)
        @variable(sp, u_sell >= 0)
        z = add_ddu_matrices(sp, [Φ(0, ρ), Φ(1, ρ)]; M = 1e4)
        if t == 1
            c_q = @constraint(sp, x.out <= x.in + 5)
            # SDDP.parameterize(ω -> set_normalized_rhs(c_q, ω), sp, Ω_q)
            @stageobjective(sp, -10 * z[2])
        else
            @constraint(sp, u_sell <= x.in)
            @constraint(sp, x.out == x.in - u_sell)
            set_upper_bound(u_sell, 10)
            # SDDP.parameterize(ω -> set_upper_bound(u_sell, ω), sp, Ω_d)
            @stageobjective(sp, 2u_sell)
        end
        return
    end
end

model = build_and_cheese_model(1234)

SDDP.train(
    model;
    forward_pass = DecisionDependentForwardPass(),
    duality_handler = # SDDP.BanditDuality(
        # SDDP.ContinuousConicDuality(),
        # SDDP.StrengthenedConicDuality(),
        SDDP.LagrangianDuality(),
    # ),
    cut_type = SDDP.MULTI_CUT,
    log_every_iteration = true,
    iteration_limit = 100,
)
ret = solve_decision_dependent_trajectory(
    model,
    model.initial_root_state,
    [:x, :u_sell];
    explore = false,
)


# # Lagrangian
# model = build_and_train_model(; Ω = Ω, Φ = Φ, z = [0, 1], ρ = 0.9)
# SDDP.train(model; forward_pass = SDDP.LagrangianDuality())

# # Lagrangian
# model = build_and_train_model(; Ω = Ω, Φ = Φ, z = [0, 1], ρ = 0.9)
# SDDP.train(model; forward_pass = SDDP.LagrangianDuality())


# model = build_and_train_model(;
#     Ω = Ω,
#     Φ = Φ,
#     z = [0, 1],
#     ρ = 0.9,
#     forward_pass =
# )

# ret = solve_decision_dependent_trajectory(
#     model,
#     model.initial_root_state,
#     [:x, :u_buy, :u_sell],
# )


# function build_and_train_model(;
#     Ω::Vector,
#     Φ::Function,
#     z::Union{Number,Vector},
#     ρ::Float64;,
# )
#     N = length(Ω)
#     graph = SDDP.Graph(0)
#     SDDP.add_node.((graph,), 1:N)
#     Φ̅ = Φ(sum(z) / length(z), ρ)
#     for j in 1:N
#         SDDP.add_edge(graph, 0 => j, Φ̅[1, j+1])
#         for i in 1:N
#             SDDP.add_edge(graph, i => j, Φ̅[i+1, j+1])
#         end
#     end
#     return SDDP.PolicyGraph(
#         graph;
#         sense = :Max,
#         upper_bound = 5 * maximum(maximum.(Ω)) / (1 - ρ),
#         optimizer = HiGHS.Optimizer,
#     ) do sp, t
#         @variable(sp, x >= 0, SDDP.State, initial_value = 0)
#         @variable(sp, u_buy >= 0)
#         @variable(sp, u_sell >= 0)
#         @constraint(sp, u_sell <= x.in)
#         @constraint(sp, x.out == x.in - u_sell + u_buy)
#         SDDP.parameterize(ω -> set_upper_bound(u_sell, ω), sp, Ω[t])
#         @stageobjective(sp, 5u_sell - 2u_buy - 0.1x.out)
#         if z isa Vector
#             _ = add_ddu_matrices(sp, Φ.(z, ρ); M = 1e-4)
#         end
#         return
#     end
# end

# function Φ(z, ρ)
#     return [
#         0.0 0.5 0.5
#         0.0 ρ * (0.5 - 0.3z) ρ * (0.5 + 0.3z)
#         0.0 ρ * (0.5 - 0.2z) ρ * (0.5 + 0.2z)
#     ]
# end

# D = [
#     Distributions.TriangularDist(150.0, 250.0, 180.0),
#     Distributions.TriangularDist(150.0, 250.0, 220.0),
# ]
# Ω = [round.(Int, sort!(rand(d, 30))) for d in D]

# # Convex relaxation
# model = build_and_train_model(; Ω = Ω, Φ = Φ, z = [0, 1], ρ = 0.9)
# SDDP.train(model; forward_pass = DecisionDependentForwardPass())

# # Lagrangian
# model = build_and_train_model(; Ω = Ω, Φ = Φ, z = [0, 1], ρ = 0.9)
# SDDP.train(model; forward_pass = SDDP.LagrangianDuality())

# # Lagrangian
# model = build_and_train_model(; Ω = Ω, Φ = Φ, z = [0, 1], ρ = 0.9)
# SDDP.train(model; forward_pass = SDDP.LagrangianDuality())


# model = build_and_train_model(;
#     Ω = Ω,
#     Φ = Φ,
#     z = [0, 1],
#     ρ = 0.9,
#     forward_pass =
# )

# ret = solve_decision_dependent_trajectory(
#     model,
#     model.initial_root_state,
#     [:x, :u_buy, :u_sell],
# )
