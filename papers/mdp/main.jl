#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP
import Gurobi
import Plots
import Random
import Statistics
import StatsPlots

###
### Section 7.1
###

function _add_ddu_constraints(model::SDDP.PolicyGraph{Int}, i::Int)
    node = model[i]
    if get(node.ext, :_ddu_is_set, false)
        return
    end
    nominal_P = [
        child.probability * noise.probability for child in node.children for
        noise in model[child.term].noise_terms
    ]
    push!(node.bellman_function.risk_set_cuts, nominal_P)
    N = length(nominal_P)
    SDDP._add_locals_if_necessary(node, node.bellman_function, N)
    θʲʷ = VariableRef[node.bellman_function.local_thetas[i].theta for i in 1:N]
    Θ = node.bellman_function.global_theta.theta
    ddu = node.subproblem.ext[:__ddu__]
    for (d, y_d) in enumerate(ddu.y)
        P_d = Float64[
            ddu.matrices[d][i+1, child.term] * noise.probability for
            child in node.children for
            noise in model[child.term].noise_terms
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
    depth::Union{Nothing,Int} = nothing,
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
        if depth === nothing &&
           rand() <=
           1 - sum(child.probability for child in node.children; init = 0)
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

function create_figure_10()
    Φ(ρ, z) = [1 0; ρ*(1-z) z; ρ 0]
    ρ = 0.9
    graph = SDDP.Graph(0)
    SDDP.add_node.((graph,), 1:2)
    Φ̅ = Φ(ρ, 0.5)
    for i in 1:3, j in 1:2
        Φ̅[i, j] > 0 && SDDP.add_edge(graph, (i - 1) => j, Φ̅[i, j])
    end
    model = SDDP.PolicyGraph(
        graph;
        sense = :Max,
        optimizer = Gurobi.Optimizer,
        upper_bound = 7 / (1 - ρ),
    ) do sp, node
        @variable(sp, x >= 0, SDDP.State, initial_value = 0)
        @variable(sp, u_sell >= 0)
        sp[:z] = z = add_ddu_matrices(sp, [Φ(ρ, 0), Φ(ρ, 1)]; M = 1e3)
        @constraint(sp, con_balance, x.out == x.in - u_sell + 0.0)
        if node == 1  # farm
            fix.(u_sell, 0; force = true)
            @stageobjective(sp, -3 * z[2])
            SDDP.parameterize(sp, [0, 2, 4, 6, 8]) do ω
                return set_normalized_rhs(con_balance, ω)
            end
        else         # market
            @stageobjective(sp, 1 * u_sell)
            @constraint(sp, u_sell <= x.in)
            SDDP.parameterize(ω -> set_upper_bound(u_sell, ω), sp, [5, 10])
        end
    end
    Random.seed!(12345)
    SDDP.train(
        model;
        duality_handler = SDDP.LagrangianDuality(),
        cut_type = SDDP.MULTI_CUT,
        forward_pass = DecisionDependentForwardPass(),
        iteration_limit = 100,
        log_every_iteration = true,
        cut_deletion_minimum = 100,
    )
    Random.seed!(5678)
    ret = solve_decision_dependent_trajectory(
        model,
        model.initial_root_state,
        [:x, :u_sell, :z];
        explore = false,
        depth = 50,
    )
    stock_plot = Plots.plot(
        map(d -> d[:x].out, ret.simulation);
        ylabel = "Quantity in stock (\$x^\\prime\$)\n",
        ylims = (0, maximum(d -> d[:x].out, ret.simulation) + 1),
        color = :slategray,
        legend = false,
        linewidth = 3,
    )
    Plots.scatter!(
        stock_plot,
        [
            (i, data[:x].out) for (i, data) in enumerate(ret.simulation) if
            data[:node_index] == 1 && data[:z][2] > 0.5
        ],
        color = "#43a047",
    )
    Plots.scatter!(
        stock_plot,
        [
            (i, data[:x].out) for (i, data) in enumerate(ret.simulation) if
            data[:node_index] == 1 && data[:z][2] < 0.5
        ],
        marker = :x,
        markerstrokewidth = 3,
        color = "#e53935",
    )
    plt = Plots.plot(
        stock_plot,
        Plots.plot(
            map(d -> d[:u_sell], ret.simulation);
            ylabel = "Sales decision (\$u_{sell}\$)",
            seriestype = :steppre,
            linewidth = 3,
            color = :slategray,
            xlabel = "Simulation step",
        ),
        xlims = (0, length(ret.simulation) + 1),
        legend = false,
        layout = (2, 1),
        dpi = 400,
    )
    Plots.savefig("cheese_producer2.pdf")
    return model, plt
end

function _solve_finite_cheese_producer(T::Int; create_plot::Bool = false)
    function Φ(z)
        a = zeros(2T + 1, 2T)
        a[1, 1] = 1.0
        for t in 1:(T-1)
            a[2t, 2t] = z
            a[2t, 2t+1] = 1.0 - z
            a[2t+1, 2t+1] = 1.0
        end
        a[2T, 2T] = z
        return a
    end
    graph = SDDP.Graph(0)
    SDDP.add_node.((graph,), 1:2T)
    Φ̅ = Φ(0.5)
    for i in 1:size(Φ̅, 1), j in 1:size(Φ̅, 2)
        Φ̅[i, j] > 0 && SDDP.add_edge(graph, (i - 1) => j, Φ̅[i, j])
    end
    model = SDDP.PolicyGraph(
        graph;
        sense = :Max,
        optimizer = Gurobi.Optimizer,
        upper_bound = 7 * T,
    ) do sp, node
        @variable(sp, x >= 0, SDDP.State, initial_value = 0)
        @variable(sp, u_sell >= 0)
        sp[:z] = z = add_ddu_matrices(sp, [Φ(0), Φ(1)]; M = 200)
        @constraint(sp, con_balance, x.out == x.in - u_sell + 0.0)
        if isodd(node)  # farm
            fix.(u_sell, 0; force = true)
            @stageobjective(sp, -3 * z[2])
            SDDP.parameterize(sp, [0, 2, 4, 6, 8]) do ω
                return set_normalized_rhs(con_balance, ω)
            end
        else         # market
            @stageobjective(sp, 1 * u_sell)
            @constraint(sp, u_sell <= x.in)
            SDDP.parameterize(ω -> set_upper_bound(u_sell, ω), sp, [5, 10])
        end
    end
    Random.seed!(123456)
    SDDP.train(
        model;
        duality_handler = SDDP.LagrangianDuality(),
        cut_type = SDDP.MULTI_CUT,
        forward_pass = DecisionDependentForwardPass(),
        iteration_limit = 100,
        log_every_iteration = true,
        cut_deletion_minimum = 500,
        stopping_rules = [SDDP.SimulationStoppingRule()],
    )
    Random.seed!(56789)
    if create_plot
        ret = solve_decision_dependent_trajectory(
            model,
            model.initial_root_state,
            [:x, :u_sell, :z];
            explore = false,
        )
        stock_plot = Plots.plot(
            map(d -> d[:x].out, ret.simulation);
            ylabel = "Quantity in stock (\$x^\\prime\$)\n",
            ylims = (0, maximum(d -> d[:x].out, ret.simulation) + 1),
            color = :slategray,
            legend = false,
            linewidth = 3,
        )
        Plots.scatter!(
            stock_plot,
            [
                (i, data[:x].out) for (i, data) in enumerate(ret.simulation) if
                isodd(data[:node_index]) && data[:z][2] > 0.5
            ],
            color = "#43a047",
        )
        Plots.scatter!(
            stock_plot,
            [
                (i, data[:x].out) for (i, data) in enumerate(ret.simulation) if
                isodd(data[:node_index]) && data[:z][2] < 0.5
            ],
            marker = :x,
            markerstrokewidth = 3,
            color = "#e53935",
        )
        plt = Plots.plot(
            stock_plot,
            Plots.plot(
                map(d -> d[:u_sell], ret.simulation);
                ylabel = "Sales decision (\$u_{sell}\$)",
                seriestype = :steppre,
                linewidth = 3,
                color = :slategray,
                xlabel = "Simulation step",
            ),
            xlims = (0, length(ret.simulation) + 1),
            legend = false,
            layout = (2, 1),
            dpi = 400,
        )
        Plots.savefig("cheese_producer_$T.pdf")
    end
    simulations = map(1:1_000) do -
        ret = solve_decision_dependent_trajectory(
            model,
            model.initial_root_state,
            Symbol[];
            explore = false,
        )
        return ret.cumulative_value
    end
    upper_bound = SDDP.calculate_bound(model)
    return upper_bound, simulations
end

function create_figure_11()
    data = Dict()
    for t in 2:2:12
        data[t] = _solve_finite_cheese_producer(t)
    end
    x = sort(collect(keys(data)))
    ub = [data[xi][1] for xi in x]
    μ = [data[xi][2] for xi in x]
    box_y = reduce(vcat, μ)
    box_x = reduce(vcat, [fill(x[i], length(μ[i])) for i in 1:length(x)])
    StatsPlots.violin(
        box_x,
        box_y;
        xlims = (1, 13),
        xticks = (2:2:12),
        xlabel = "Time horizon",
        ylabel = "Objective value",
        label = false,
        color = :grey,
        alpha = 0.5,
    )
    Plots.scatter!(
        x,
        ub;
        label = "Upper bound",
        color = :black,
    )
    Plots.scatter!(
        x,
        Statistics.mean.(μ);
        label = "Sample mean",
        marker = :o,
        color = :white,
    )
    Plots.savefig("cheese_producer_violin.pdf")
    return
end

###
### Section 7.2
###

mutable struct StoppingForwardPass <: SDDP.AbstractForwardPass
    name::Symbol
end

function SDDP.forward_pass(
    model::SDDP.PolicyGraph,
    options::SDDP.Options,
    f_pass::StoppingForwardPass,
)
    pass = SDDP.forward_pass(model, options, SDDP.DefaultForwardPass())
    index = findfirst(s -> s[f_pass.name] < 0.5, pass.sampled_states)
    n = 3 + something(index, length(pass.sampled_states))
    subset(x, n) = x[1:min(length(x), n)]
    return (
        scenario_path = subset(pass.scenario_path, n),
        sampled_states = subset(pass.sampled_states, n),
        objective_states = pass.objective_states,
        belief_states = subset(pass.belief_states, n),
        cumulative_value = pass.cumulative_value,
    )
end

function _solve_tiger_problem(ε::Float64; create_plot::Bool = false)
    ρ = 0.95
    graph = SDDP.Graph(
        :R,
        [:l, :r],
        [(:R => :l, 0.5), (:R => :r, 0.5), (:l => :l, ρ), (:r => :r, ρ)],
    )
    if ε < 0.5
        SDDP.add_ambiguity_set(graph, [:l, :r], 1e3)
    end
    model = SDDP.PolicyGraph(
        graph;
        sense = :Min,
        lower_bound = -10.0,
        optimizer = Gurobi.Optimizer,
    ) do sp, node
        # s: stay, l: open left, r: open right
        @variable(sp, x_s, Bin, SDDP.State, initial_value = 1)
        @variable(sp, x_l, Bin, SDDP.State, initial_value = 0)
        @variable(sp, x_r, Bin, SDDP.State, initial_value = 0)
        @constraint(sp, x_s.out + x_l.out + x_r.out <= 1 - x_l.in - x_r.in)
        @constraint(sp, x_s.out + x_l.out + x_r.out == x_s.in)
        if node == :l
            @stageobjective(sp, 100 * x_l.in - 10 * x_r.in + x_s.in)
            SDDP.parameterize(sp, [:left, :right], [0.5 + ε, 0.5 - ε]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        elseif node == :r
            @stageobjective(sp, -10 * x_l.in + 100 * x_r.in + x_s.in)
            SDDP.parameterize(sp, [:left, :right], [0.5 - ε, 0.5 + ε]) do ω
                # println("I heard the tiger on the $ω side.")
            end
        end
    end
    Random.seed!(12345)
    SDDP.train(
        model;
        iteration_limit = 100,
        log_every_iteration = true,
        cut_deletion_minimum = 1_000,
        duality_handler = SDDP.LagrangianDuality(),
        forward_pass = StoppingForwardPass(:x_s),
    )
    lower_bound = SDDP.calculate_bound(model)
    Random.seed!(4567)
    sampling_scheme =
        SDDP.InSampleMonteCarlo(max_depth = 50, terminate_on_dummy_leaf = false)
    simulations =
        SDDP.simulate(model, 100, [:x_s, :x_l, :x_r]; sampling_scheme)
    objectives = map(simulations) do simulation
        return sum(
            ρ^(t - 1) * d[:stage_objective] for (t, d) in enumerate(simulation)
        )
    end
    if !create_plot
        μ, σ = SDDP.confidence_interval(objectives)
        println("lower_bound = $lower_bound")
        println("upper_bound = $μ ± $σ")
        return lower_bound, objectives, model
    end
    simulations = simulations[1:100]
    belief_plot = Plots.plot(;
        xlabel = "Time step",
        ylabel = "Belief(Left)",
        legend = false,
        ymajorgrid = true,
    )
    plot = Plots.plot(;
        xlabel = "Time step",
        ylabel = "# hear left - hear right",
        legend = false,
        ylims = (-4, 4),
        ymajorgrid = true,
    )
    for simulation in simulations
        b = Float64[0.5]
        y = Int[0]
        for d in simulation
            push!(b, d[:belief][:l])
            if d[:noise_term] == :left
                push!(y, y[end] + 1)
            else
                push!(y, y[end] - 1)
            end
            if d[:x_l].out > 0.5 || d[:x_r].out > 0.5
                break
            end
        end
        Plots.plot!(
            belief_plot,
            0:length(b)-1,
            b;
            color = :grey,
            linewidth = 3,
            alpha = 0.2,
        )
        Plots.plot!(
            plot,
            0:length(y)-1,
            y;
            color = :grey,
            linewidth = 3,
            alpha = 0.2,
        )
        function correct_door(d)
            return (d[:x_l].out > 0.5 && d[:node_index] == :r) ||
                   (d[:x_r].out > 0.5 && d[:node_index] == :l)
        end
        function incorrect_door(d)
            return (d[:x_l].out > 0.5 && d[:node_index] == :l) ||
                   (d[:x_r].out > 0.5 && d[:node_index] == :r)
        end
        if (i = findfirst(correct_door, simulation)) !== nothing
            Plots.scatter!([i], [y[i+1]], color = "#43a047", markersize = 6)
        end
        if (i = findfirst(incorrect_door, simulation)) !== nothing
            @show i
            Plots.scatter!(
                [i],
                [y[i+1]];
                marker = :x,
                markersize = 6,
                markerstrokewidth = 3,
                color = "#e53935",
            )
        end
    end
    Plots.plot(belief_plot, plot, layout = (2, 1), dpi = 400)
    Plots.savefig("tiger_problem_$ε.pdf")
    return lower_bound, objectives, model
end

create_figure_12() = _solve_tiger_problem(0.35; create_plot = true)

function create_figure_13()
    data = Dict()
    for ε in [0.2, 0.3, 0.4, 0.5]
        data[ε] = _solve_tiger_problem(ε)
    end

    x = sort(collect(keys(data)))
    lb = [data[xi][1] for xi in x]
    μ = [data[xi][2] for xi in x]
    box_y = reduce(vcat, μ)
    box_x = reduce(vcat, [fill(0.5 - x[i], length(μ[i])) for i in 1:length(x)])
    StatsPlots.violin(
        box_x,
        # Apply small perturbationn for the violin to show up when constant.
        box_y .+ 0.15 * rand(length(box_y));
        bar_width = 0.05,
        # ylims = (-10, 25),
        xlabel = "False positive rate",
        ylabel = "Objective value",
        label = false,
        color = :grey,
        alpha = 0.5,
    )
    Plots.scatter!(
        0.5 .- x,
        lb;
        label = "Lower bound",
        color = :black,
        linewidth = 3,
    )
    Plots.scatter!(
        0.5 .- x,
        Statistics.mean.(μ);
        label = "Sample mean",
        marker = :o,
        color = :white,
    )
    Plots.savefig("tiger_problem_violin.pdf")
    return
end

# create_figure_10()
# create_figure_11()
# create_figure_12()
# create_figure_13()
