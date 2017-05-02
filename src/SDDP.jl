#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

module SDDP

using JuMP, Distributions

export SDDPModel,
    # inputs
    @state, @states,
    @scenario, @scenarios, setscenarioprobability!,
    stageobjective!,
    # cut oracles
    DefaultCutOracle,
    # risk measures
    Expectation,
    solve

include("storage.jl")
include("typedefinitions.jl")
include("utilities.jl")
include("riskmeasures.jl")
include("states.jl")
include("scenarios.jl")
include("cutoracles.jl")
include("valuefunctions.jl")
include("default_value_function.jl")
include("stageobjectives.jl")
include("cuttingpasses.jl")
include("MIT_licensedcode.jl")
include("print.jl")

include("pro/dematos_cutselection.jl")
include("pro/avar_riskaversion.jl")
include("pro/rib_valuefunction.jl")

struct UnsetSolver <: JuMP.MathProgBase.AbstractMathProgSolver end

function SDDPModel(build!::Function;
    sense           = :Min,
    stages          = 1,
    objective_bound = -1e6,

    # markov_states        = 1,
    # initial_markov_state = -1,
    markov_transition    = [ones(Float64, (1,1)) for t in 1:stages],

    scenario_probability = Float64[],

    risk_measure    = Expectation(),

    cut_oracle      = DefaultCutOracle(),

    solver          = UnsetSolver(),

    value_function = DefaultValueFunction(),
    )

    # check number of arguments to SDDPModel() do [args...] ...  model def ... end
    num_args = n_args(build!)
    if num_args < 2
        error("""Too few arguments in
            SDDPModel() do args...
            end""")
    elseif num_args > 4
        error("""Too many arguments
            SDDPModel() do args...
            end""")
    end

    # New SDDPModel
    m = newSDDPModel(value_function, build!, solver)

    for t in 1:stages
        # todo: transition probabilities that vary by stage
        stage = Stage(getel(Array{Float64, 2}, markov_transition, t))
        # check that
        if num_args == 2 && size(markov_transition[t], 2) > 1
            error("""Because you specified a scenario tree in the SDDPModel constructor, you need to use the

                SDDPModel() do sp, stage, markov_state
                    ... model definition ...
                end

            syntax.""")
        end
        for i in 1:size(markov_transition[t], 2)
            mod = Subproblem(
                finalstage     = (t == stages),
                stage          = t,
                markov_state   = i,
                sense          = optimisationsense(sense),
                bound          = float(objective_bound),
                risk_measure   = getel(AbstractRiskMeasure, risk_measure, t, i),
                cut_oracle     = deepcopy(getel(AbstractCutOracle, cut_oracle, t, i)),
                value_function = value_function
            )
            setsolver(mod, solver)
            # dispatch to correct function
            # maybe we should do this with tuples
            if num_args == 3
                build!(mod, t, i)
            else
                build!(mod, t)
            end
            # Uniform scenario probability for now
            scenario_prob = getel(Vector, scenario_probability, t, i)
            if length(scenario_prob) != 0 && length(scenario_prob) != length(ext(mod).scenarios)
                error("Invalid number of scenarios.")
            end
            if length(scenario_prob) == 0
                for i in 1:length(ext(mod).scenarios)
                    push!(ext(mod).scenarioprobability, 1 / length(ext(mod).scenarios))
                end
            else
                if abs(sum(scenario_prob) - 1) > 1e-6 # check probability
                    error("You must specify a discrete probability distribution that sums to 1.0")
                end
                for i in 1:length(ext(mod).scenarios)
                    push!(ext(mod).scenarioprobability, scenario_prob[i])
                end
            end
            push!(stage.subproblems, mod)
        end
        push!(m.stages, stage)
    end
    m
end

contains(x, l, u) = x >= l && x <= u

applicable(iteration, frequency) = frequency > 0 && mod(iteration, frequency) == 0

function solve_serial(m::SDDPModel, settings::Settings=Settings())
    status = :solving
    time_simulating, time_cutting = 0.0, 0.0
    objectives = CachedVector(Float64)
    nsimulations, iteration, keep_iterating = 0, 1, true
    start_time = time()
    while keep_iterating
        (objective_bound, time_backwards, simulation_objective, time_forwards) = iteration!(m, settings)
        time_cutting += time_backwards + time_forwards
        lower, upper = simulation_objective, simulation_objective

        # cut selection
        if applicable(iteration, settings.cut_selection_frequency)
            settings.print_level > 1 && info("Running Cut Selection")
            for (t, stage) in enumerate(stages(m))
                t == length(stages(m)) && continue
                for sp in subproblems(stage)
                    rebuildsubproblem!(m, sp)
                end
            end
        end

        # simulate policy
        if applicable(iteration, settings.simulation_frequency)
            settings.print_level > 1 && info("Running Monte-Carlo Simulation")
            t = time()
            reset!(objectives)
            simidx = 1
            for i in 1:settings.simulation_steps[end]
                push!(objectives, forwardpass!(m, settings))
                nsimulations += 1
                if i == settings.simulation_steps[simidx]
                    (lower, upper) = confidenceinterval(objectives, settings.simulation_confidence)
                    if contains(objective_bound, lower, upper)
                        if settings.simulation_terminate && simidx == length(settings.simulation_steps)
                            # max simulations
                            status = :converged
                            keep_iterating = false
                        end
                    else
                        break
                    end
                    simidx += 1
                end
            end
            time_simulating += time() - t
        end

        total_time = time() - start_time
        push!(m.log, SolutionLog(iteration, objective_bound, lower, upper, time_cutting, nsimulations, time_simulating, total_time))

        settings.print_level > 0 && print(STDOUT, m.log[end], mod(iteration, settings.simulation_frequency) != 0)
        settings.log_file != "" && print(settings.log_file, m.log[end], mod(iteration, settings.simulation_frequency) != 0)

        iteration += 1
        if iteration > settings.max_iterations
            status = :max_iterations
            keep_iterating = false
        end

    end
end

function solve(m::SDDPModel;
        max_iterations::Int       = 10,
        simulation_frequency::Int = 5,
        simulation_min::Int       = 10,
        simulation_max::Int       = 20,
        simulation_step::Int      = 10,
        simulation_confidence::Float64 = 0.95,
        simulation_terminate::Bool = false,
        cut_selection_frequency::Int = 0,
        print_level::Int          = 4,
        log_file::String          = ""
    )
    settings = Settings(
        max_iterations,
        simulation_frequency,
        collect(simulation_min:simulation_step:simulation_max),
        simulation_confidence,
        simulation_terminate,
        cut_selection_frequency,
        print_level,
        log_file
    )

    print_level > 0 && printheader(STDOUT, m)
    log_file != "" && printheader(log_file, m)

    solve_serial(m, settings)

end

end
