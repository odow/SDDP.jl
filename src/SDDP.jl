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

function SDDPModel(build!::Function;
    sense           = :Min,
    stages          = 1,
    objective_bound = -1e6,
    markov_states   = 1,
    initial_markov_probability = [1.0],
    transition      = Array{Float64}(0,0),
    risk_measure    = Expectation(),
    cut_oracle      = DefaultCutOracle(),
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
                stage        = t,
                markov_state = i,
                sense        = optimisationsense(sense),
                bound        = float(objective_bound),
                risk_measure = getel(AbstractRiskMeasure, risk_measure, t, i),
                cut_oracle   = getel(AbstractCutOracle, cut_oracle, t, i),
                value_function = DefaultValueFunction
            )
            # dispatch to correct function
            # maybe we should do this with tuples
            if num_args == 3
                build!(mod, t, i)
            else
                build!(mod, t)
            end
            # Uniform scenario probability for now
            for i in 1:length(ext(mod).scenarios)
                push!(ext(mod).scenarioprobability, 1 / length(ext(mod).scenarios))
            end
            push!(stage.subproblems, mod)
        end
        push!(m.stages, stage)
    end
    m
end

function solve_serial(m::SDDPModel, settings::Settings=Settings())
    for itr in 1:settings.max_iterations
        iteration!(m, settings)
        # simulate policy
        if mod(itr, settings.simulation_frequency) == 0
            #
        end


    end
end

end
