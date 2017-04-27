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

using JuMP

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
include("stageobjectives.jl")
include("cuttingpasses.jl")
include("MIT_licensedcode.jl")
include("print.jl")

struct UnsetSolver <: JuMP.MathProgBase.AbstractMathProgSolver end

function SDDPModel(build!::Function;
    sense           = :Min,
    stages          = 1,
    objective_bound = -1e6,

    markov_states   = 1,
    initial_markov_probability = ones(Float64, 1),
    transition      = ones(Float64, (1,1)),

    scenario_probability = Float64[],

    risk_measure    = Expectation(),

    cut_oracle      = DefaultCutOracle(),

    solver          = UnsetSolver(),

    # value_function = DefaultValueFunction{cut_oracle},

    kwargs...)

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
    m = SDDPModel()

    for t in 1:stages
        # todo: transition probabilities that vary by stage
        stage = Stage(transition)
        # check that
        if num_args == 2 && getel(Int, markov_states, t) != 1
            error("""Because you specified a scenario tree in the SDDPModel constructor, you need to use the

                SDDPModel() do sp, stage, markov_state
                    ... model definition ...
                end

            syntax.""")
        end
        for i in 1:getel(Int, markov_states, t)
            mod = Subproblem(
                finalstage   = (t == stages),
                stage        = t,
                markov_state = i,
                sense        = optimisationsense(sense),
                bound        = float(objective_bound),
                risk_measure = getel(AbstractRiskMeasure, risk_measure, t, i),
                cut_oracle   = getel(AbstractCutOracle, cut_oracle, t, i),
                value_function = DefaultValueFunction
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

function solve_serial(m::SDDPModel, settings::Settings=Settings())
    time_simulating = 0.0
    time_cutting    = 0.0
    for itr in 1:settings.max_iterations
        (objective_bound, time_backwards, simulation_objective, time_forwards) = iteration!(m, settings)
        time_cutting += time_backwards + time_forwards
        # simulate policy
        nsimulations = 0
        if mod(itr, settings.simulation_frequency) == 0
            t = time()
            #
            time_simulating += time() - t
        end

        push!(m.log, SolutionLog(itr, objective_bound, simulation_objective, simulation_objective, time_cutting, nsimulations, time_simulating))

        settings.print_level > 0 && print(STDOUT, m.log[end])
        settings.log_file != "" && print(settings.log_file, m.log[end])

    end
end

function solve(m::SDDPModel;
        max_iterations::Int       = 10,
        simulation_frequency::Int = 1,
        print_level::Int          = 4,
        log_file::String          = ""
    )
    settings = Settings(max_iterations, simulation_frequency, print_level, log_file)

    print_level > 0 && printheader(STDOUT)
    log_file != "" && printheader(log_file)

    solve_serial(m, settings)

end




end
