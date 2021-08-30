#  Copyright 2017-21, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

include(joinpath(@__DIR__, "BiObjectiveSDDP.jl"))
using .BiObjectiveSDDP

using SDDP
import Gurobi
import Statistics

const OBJ_1_SCALING = 0.01
const OBJ_2_SCALING = 0.1
include("brazilian_data.jl")

function create_model(weight = nothing)
    env = Gurobi.Env()
    model = SDDP.LinearPolicyGraph(
        stages = 12,
        lower_bound = 0.0,
        optimizer = () -> Gurobi.Optimizer(env),
    ) do sp, t
        set_silent(sp)
        month = t % 12 == 0 ? 12 : t % 12  # Year to month conversion.
        @variables(
            sp,
            begin
                0 <= storedEnergy[i = 1:4] <= storedEnergy_ub[i],
                (SDDP.State, initial_value = storedEnergy_initial[i])
                0 <= spillEnergy[i = 1:4]
                0 <= hydroGeneration[i = 1:4] <= hydro_ub[i]
                thermal_lb[i][j] <=
                thermal[i = 1:4, j = 1:N_THERMAL[i]] <=
                thermal_ub[i][j]
                0 <= exchange[i = 1:5, j = 1:5] <= exchange_ub[i][j]
                0 <=
                deficit[i = 1:4, j = 1:4] <=
                demand[month][i] * deficit_ub[j]
                inflow[i = 1:4]
            end
        )
        @constraints(
            sp,
            begin
                [i = 1:4],
                sum(deficit[i, :]) +
                hydroGeneration[i] +
                sum(thermal[i, j] for j in 1:N_THERMAL[i]) +
                sum(exchange[:, i]) - sum(exchange[i, :]) ==
                demand[month][i]
                [i = 1:4],
                storedEnergy[i].out + spillEnergy[i] + hydroGeneration[i] -
                storedEnergy[i].in == inflow[i]
                sum(exchange[:, 5]) == sum(exchange[5, :])
            end
        )
        Ω = if t == 1
            [inflow_initial]
        else
            r = (t - 1) % 12 == 0 ? 12 : (t - 1) % 12
            [
                [scenarios[i][r][ω] for i in 1:4] for
                ω in 1:length(scenarios[1][r])
            ]
        end
        @expressions(
            sp,
            begin
                objective_1,
                OBJ_1_SCALING *
                sum(deficit_obj[i] * sum(deficit[i, :]) for i in 1:4)
                objective_2,
                OBJ_2_SCALING * sum(
                    thermal_obj[i][j] * thermal[i, j] for i in 1:4 for
                    j in 1:N_THERMAL[i]
                )
            end
        )
        if weight === nothing
            SDDP.initialize_biobjective_subproblem(sp)
        end
        SDDP.parameterize(sp, Ω) do ω
            JuMP.fix.(inflow, ω)
            if weight === nothing
                SDDP.set_biobjective_functions(sp, objective_1, objective_2)
            else
                @stageobjective(
                    sp,
                    weight * objective_1 + (1 - weight) * objective_2,
                )
            end
            return
        end
    end
    return model
end

function _simulate_policy(model, keys)
    simulations = Dict()
    for λ in keys
        BiObjectiveSDDP.set_scalarizing_weight(model, λ)
        simulations[λ] = SDDP.simulate(
            model,
            1000,
            [:storedEnergy, :objective_1, :objective_2],
        )
    end
    return simulations
end

function _extract_objectives(simulation)
    obj_1 = [sum(s[:objective_1] for s in sim) for sim in simulation]
    obj_2 = [sum(s[:objective_2] for s in sim) for sim in simulation]
    return obj_1, obj_2
end

function _save_simulations_to_dat(simulations, simulation_weights)
    A = Matrix{Float64}(
        undef,
        length(simulations[first(simulation_weights)]),
        2 * length(simulation_weights),
    )
    for (i, weight) in enumerate(simulation_weights)
        obj_1, obj_2 = _extract_objectives(simulations[weight])
        A[:, 2*i-1] .= obj_1
        A[:, 2*i] .= obj_2
    end
    open("simulations.dat", "w") do io
        for i in 1:size(A, 1)
            println(io, join(A[i, :], "  "))
        end
    end
    open("new_data.dat", "w") do io
        for (i, w) in enumerate([0.1, 0.7, 0.9])
            s = w .* A[:, 2i-1] + (1 - w) .* A[:, 2i]
            μ = Statistics.mean(s)
            Q = Statistics.quantile(s, [0.1, 0.9])
            println(io, w, " ", μ, " ", μ - Q[1], " ", Q[2] - μ)
        end
    end
    return
end

"""
    experiment_1()

Run the first experiment where we train the policy using the true biobjective
SDDP algorithm.
"""
function experiment_1()
    model = create_model()
    env = Gurobi.Env()
    lower_bound, weights, bounds = BiObjectiveSDDP.bi_objective_sddp(
        model,
        () -> Gurobi.Optimizer(env);
        # BiObjectiveSDDP kwargs ...
        bi_objective_sddp_iteration_limit = 2000,
        bi_objective_lambda_atol = 0.05,
        bi_objective_lower_bound = 0.0,
        # SDDP.jl kwargs ...
        print_level = 0,
        iteration_limit = 1,
    )
    open("bounds.dat", "w") do io
        for (w, b) in zip(weights, bounds)
            println(io, "$(w)  $(b)")
        end
    end
    simulation_weights = [0.1, 0.7, 0.9]
    simulations = _simulate_policy(model, simulation_weights);
    _save_simulations_to_dat(simulations, simulation_weights)
    return
end

"""
    BoundLimit(limit::Float64)

Terminate once the bound is better than `limit`.
"""
struct BoundLimit <: SDDP.AbstractStoppingRule
    limit::Float64
    atol::Float64
end

SDDP.stopping_rule_status(::BoundLimit) = :bound_limit

function SDDP.convergence_test(
    model::SDDP.PolicyGraph,
    log::Vector{SDDP.Log},
    rule::BoundLimit,
)
    if model.objective_sense == MOI.MIN_SENSE
        return log[end].bound >= rule.limit - rule.atol
    else
        return log[end].bound <= rule.limit + rule.atol
    end
end

"""
    experiment_2(N::Int, atol::Float64)

Run an experiment in which we time how long it takes to solve N different
policies.
"""
function experiment_2(N::Int, atol::Float64)
    # Precompilation to avoid measuring that overhead!
    _model = create_model(1.0)
    SDDP.train(_model; iteration_limit = 1, print_level = 0)
    # Now the real model
    weights = [0.0, 1.0]
    queue = [(0.0, 1.0)]
    while length(weights) < N
        (a, b) = popfirst!(queue)
        c = (a + b) / 2
        push!(weights, c)
        push!(queue, (a, c))
        push!(queue, (c, b))
    end
    start_time = time()
    for weight in weights
        model = create_model(weight)
        SDDP.train(
            model;
            log_file = "experiment_2_$(weight).txt",
            stopping_rules = [SDDP.BoundStalling(10, atol)],
            # Turn of cut selection for this experiment. We don't have it for
            # the interpolation stuff.
            cut_deletion_minimum = 10_000,

        )
        bound = SDDP.calculate_bound(model)
        open("experiment_2.dat", "a") do io
            println(io, weight, ", ", bound, ", ", time() - start_time)
        end
    end
    return
end

"""
    experiment_3(atol::Float64)

Run an experiment in which we time how long it takes to solve the problems from
experiment_2 using the saddle cuts.
"""
function experiment_3(atol::Float64)
    # Precompilation to avoid measuring that overhead!
    _model = create_model()
    SDDP.train_biobjective(
        _model;
        solution_limit = 1,
        iteration_limit = 1,
        print_level = 0,
    )
    # Now the real model
    limit_pairs = Pair{Float64,BoundLimit}[]
    open("experiment_2.dat", "r") do io
        for line in readlines(io)
            items = parse.(Float64, String.(split(line, ",")))
            push!(limit_pairs, items[1] => BoundLimit(items[2], atol))
        end
    end
    limit_dict = Dict(limit_pairs)
    model = create_model()
    solutions = SDDP.train_biobjective(
        model;
        solution_limit = length(limit_pairs),
        include_timing = true,
        print_level = 1,
        log_file_prefix = "experiment_3",
        stopping_rules = weight -> [limit_dict[weight]],
    )
    open("experiment_3.dat", "w") do io
        for (weight, _) in limit_pairs
            bound, time = solutions[weight]
            println(io, weight, ", ", bound, ", ", time)
        end
    end
    return
end

function arg(T, key)
    i = findfirst(isequal(key), ARGS)
    return i === nothing ? nothing : parse(T, ARGS[i+1])
end

function help()
    println("""julia brazilian_example.jl --experiment={1,2,3} -n 9 -atol 100

    ## Arguments

     * --experiment :: choose which experiment to Run

    ### Experiment 2

     * -n    :: Choose how many weights to run
     * -atol :: The tolerance used by the `BoundStalling` stopping rule

    ## Examples

    ```
    nohup ~/julia1.6 --project=. brazilian_example.jl --experiment=1 &
    nohup ~/julia1.6 --project=. brazilian_example.jl --experiment=2 -n 9 -atol 1e2 &
    nohup ~/julia1.6 --project=. brazilian_example.jl --experiment=3 &
    ```
    """)
end

function main()
    if findfirst(isequal("--experiment=1"), ARGS) !== nothing
        experiment_1()
    elseif findfirst(isequal("--experiment=2"), ARGS) !== nothing
        experiment_2(
            something(arg(Int, "-n"), 9),
            something(arg(Float64, "-atol"), 1e2),
        )
    elseif findfirst(isequal("--experiment=3"), ARGS) !== nothing
        experiment_3(something(arg(Float64, "-atol"), 1e2))
    else
        help()
    end
    return
end

if length(ARGS) > 0
    main()
end
