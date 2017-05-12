#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
    MonteCarloSimulation, BoundConvergence,
    Serial,
    solve, simulate,
    getbound

include("typedefinitions.jl")
include("utilities.jl")
include("riskmeasures.jl")
include("states.jl")
include("scenarios.jl")
include("cutoracles.jl")
include("valuefunctions.jl")
include("simulate.jl")
include("MIT_licensedcode.jl")
include("print.jl")

include("pro/dematos_cutselection.jl")
include("pro/avar_riskaversion.jl")
include("pro/rib_valuefunction.jl")
include("pro/solve_asyncronous.jl")
include("pro/visualiser/visualise.jl")

struct UnsetSolver <: JuMP.MathProgBase.AbstractMathProgSolver end

function SDDPModel(build!::Function;
    sense                = :Min,
    stages               = 1,
    objective_bound      = -1e6,
    markov_transition    = [ones(Float64, (1,1)) for t in 1:stages],
    scenario_probability = Float64[],
    risk_measure         = Expectation(),
    cut_oracle           = DefaultCutOracle(),
    solver               = UnsetSolver(),
    value_function       = DefaultValueFunction(cut_oracle),
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
    m = newSDDPModel(sense, value_function, build!, solver)

    for t in 1:stages
        # todo: transition probabilities that vary by stage
        stage = Stage(t, getel(Array{Float64, 2}, markov_transition, t))
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
                value_function = deepcopy(value_function)
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

samplesubproblem(stage::Stage, last_markov_state::Int, solutionstore::Void) = samplesubproblem(stage, last_markov_state)

function samplesubproblem(stage::Stage, last_markov_state, solutionstore::Dict{Symbol, Any})
    if length(solutionstore[:scenario]) >= stage.t
        idx = solutionstore[:markov][stage.t]
        return idx, getsubproblem(stage, idx)
    else
        return samplesubproblem(stage, last_markov_state)
    end
end

samplescenario(sp::JuMP.Model, solutionstore::Void) = samplescenario(sp)

function samplescenario(sp::JuMP.Model, solutionstore::Dict{Symbol, Any})
    if length(solutionstore[:scenario])>=ext(sp).stage
        idx = solutionstore[:scenario][ext(sp).stage]
        return idx, ext(sp).scenarios[idx]
    else
        return samplescenario(sp)
    end
end

function forwardpass!(m::SDDPModel, settings::Settings, solutionstore=nothing)
    last_markov_state = 1
    scenarioidx = 0
    obj = 0.0
    for (t, stage) in enumerate(stages(m))
        # choose markov state
        (last_markov_state, sp) = samplesubproblem(stage, last_markov_state, solutionstore)
        if t > 1
            setstates!(m, sp)
        end

        # choose and set RHS scenario
        if hasscenarios(sp)
            (scenarioidx, scenario) = samplescenario(sp, solutionstore)
            setscenario!(sp, scenario)
        end
        # solve subproblem
        solvesubproblem!(ForwardPass, m, sp)
        # store stage obj
        obj += getstageobjective(sp)
        # store state
        savestates!(stage.state, sp)
        # save solution for simulations (defaults to no-op)
        savesolution!(solutionstore, last_markov_state, scenarioidx, sp, t)
    end
    return obj
end

function backwardpass!(m::SDDPModel, settings::Settings)
    # walk backward through the stages
    for t in nstages(m):-1:2
        # solve all stage t problems
        reset!(m.storage)
        for sp in subproblems(m, t)
            setstates!(m, sp)
            solvesubproblem!(BackwardPass, m, sp)
        end
        # add appropriate cuts
        for sp in subproblems(m, t-1)
            modifyvaluefunction!(m, settings, sp)
        end
    end

    reset!(m.storage)
    for sp in subproblems(m, 1)
        solvesubproblem!(BackwardPass, m, sp)
    end
    # TODO: improve over just taking mean of first stage subproblems
    bound = mean(m.storage.objective)

    return bound
end

function iteration!(m::SDDPModel, settings::Settings)
    t = time()

    simulation_objective = forwardpass!(m, settings)
    time_forwards = time() - t

    objective_bound = backwardpass!(m, settings)
    time_backwards = time() - time_forwards - t

    # reduce memory footprint
    if settings.reduce_memory_footprint
        # For future reference
        # https://github.com/JuliaOpt/JuMP.jl/issues/969#issuecomment-282191105
        for stage in stages(m)
            for sp in subproblems(stage)
                for con in sp.linconstr
                    con.terms = 0
                end
            end
        end
    end
    return objective_bound, time_backwards, simulation_objective, time_forwards

end

function JuMP.solve(::Serial, m::SDDPModel, settings::Settings=Settings())
    status = :solving
    time_simulating, time_cutting = 0.0, 0.0
    objectives = CachedVector(Float64)
    nsimulations, iteration, keep_iterating = 0, 1, true
    start_time = time()
    while keep_iterating
        # add cuts
        (objective_bound, time_backwards, simulation_objective, time_forwards) = iteration!(m, settings)
        # update timers and bounds
        time_cutting += time_backwards + time_forwards
        lower, upper = simulation_objective, simulation_objective

        if applicable(iteration, settings.cut_selection_frequency)
            # run cut selection
            settings.print_level > 1 && info("Running Cut Selection")
            for (t, stage) in enumerate(stages(m))
                t == length(stages(m)) && continue
                for sp in subproblems(stage)
                    rebuildsubproblem!(m, sp)
                end
            end
        end

        if applicable(iteration, settings.simulation.frequency)
            # simulate policy
            settings.print_level > 1 && info("Running Monte-Carlo Simulation")
            t = time()
            simidx = 1
            # reuse store for objectives
            reset!(objectives)
            for i in 1:settings.simulation.steps[end]
                # forwardpass! returns objective
                push!(objectives, forwardpass!(m, settings))
                nsimulations += 1
                if i == settings.simulation.steps[simidx]
                    # simulation incrementation
                    (lower, upper) = confidenceinterval(objectives, settings.simulation.confidence)
                    if contains(objective_bound, lower, upper)
                        if settings.simulation.termination && simidx == length(settings.simulation.steps)
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

        addsolutionlog!(m, settings, iteration, objective_bound, lower, upper, time_cutting, nsimulations, time_simulating, total_time, !applicable(iteration, settings.simulation.frequency))

        status, keep_iterating = testconvergence(m, settings)

        if total_time > settings.time_limit
            status = :time_limit
            keep_iterating = false
        end

        iteration += 1
        if iteration > settings.max_iterations
            status = :max_iterations
            keep_iterating = false
        end

    end
    status
end

function addsolutionlog!(m, settings, iteration, objective, lower, upper, cutting_time, simulations, simulation_time, total_time, printsingle)
    push!(m.log, SolutionLog(iteration, objective, lower, upper, cutting_time, simulations, simulation_time, total_time))
    print(print, settings, m.log[end], printsingle)
end

function testconvergence(m::SDDPModel, settings::Settings)
    if settings.bound_convergence.iterations > 1 && length(m.log) >= settings.bound_convergence.iterations
        last_n = map(l->l.bound, m.log[end-settings.bound_convergence.iterations+1:end])
        if all(last_n - mean(last_n) .< settings.bound_convergence.atol) || all(abs.(last_n / mean(last_n)-1) .<    settings.bound_convergence.rtol)
            return :bound_convergence, false
        end
    end
    return :solving, true
end

function JuMP.solve(m::SDDPModel;
        max_iterations::Int       = 10,
        time_limit::Real          = 600, # seconds
        simulation = MonteCarloSimulation(
                frequency   = 0,
                min         = 10,
                step        = 10,
                max         = 20,
                confidence  = 0.95,
                termination = false
            ),
        bound_convergence = BoundConvergence(
                iterations = 0,
                rtol       = 0.0,
                atol       = 0.0
            ),
        cut_selection_frequency::Int = 0,
        print_level::Int             = 4,
        log_file::String             = "",
        solve_type::SDDPSolveType    = Serial(),
        # this reduces memory but you shouldn't use it if you want to save the
        # sddp model since it throws away some information
        reduce_memory_footprint      = false,
        cut_output_file::String      = ""
    )
    settings = Settings(
        max_iterations,
        time_limit,
        simulation,
        bound_convergence,
        cut_selection_frequency,
        print_level,
        log_file,
        reduce_memory_footprint,
        cut_output_file
    )

    if cut_output_file != ""
        # clear it
        open(cut_output_file, "w") do file
        end
    end

    print(printheader, settings, m, solve_type)
    status = solve(solve_type, m, settings)
    print(printfooter, settings, m, status)

    status
end

end
