#  Copyright 2017, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

__precompile__()

module SDDP

using JuMP, Distributions, JSON, TimerOutputs

const TIMER = TimerOutput()

const JuMPVERSION = Pkg.installed("JuMP")

export SDDPModel,
    # inputs
    @state, @states,
    @rhsnoise, @rhsnoises, setnoiseprobability!,
    stageobjective!, @stageobjective,
    # cut oracles
    DefaultCutOracle, DematosCutOracle, LevelOneCutOracle,
    MonteCarloSimulation, BoundConvergence,
    Serial, Asynchronous,
    solve, simulate,
    # @visualise,
    getbound,
    loadcuts!

include("JuMPfunctions.jl")
include("typedefinitions.jl")
include("utilities.jl")
include("risk_measures/riskmeasures.jl")
include("states.jl")
include("noises.jl")
include("cut_oracles/cutoracles.jl")
include("defaultvaluefunction.jl")
include("simulate.jl")
include("MIT_licensedcode.jl")
include("print.jl")
include("solve_asynchronous.jl")
include("visualizer/visualize.jl")
include("price_interpolation/price_interpolation.jl")
include("deprecate.jl")

struct UnsetSolver <: JuMP.MathProgBase.AbstractMathProgSolver end

"""
    SDDPModel(;kwargs...) do ...

    end

# Description

This function constructs an SDDPModel. `SDDPModel` takes the following keyword
arguments. Some are required, and some are optional.

# Required Keyword arguments

 * `stages::Int`
 The number of stages in the problem. A stage is defined as each step in time at
 which a decion can be made. Defaults to `1`.

 * `objective_bound`
 A valid bound on the initial value/cost to go. i.e. for maximisation this may
 be some large positive number, for minimisation this may be some large negative
 number. Users can pass either a single value (which bounds the cost-to-go in all
 stages), or a vector of values (one for each stage), or a vector (one element
 for each stage) of vectors of values (one value for each markov state in the stage).

 * `solver::MathProgBase.AbstractMathProgSolver`
 MathProgBase compliant solver that returns duals from a linear program. If this
 isn't specified then you must use `JuMP.setsolver(sp, solver)` in the stage
 definition.

# Optional Keyword arguments

 * `sense`
 Must be either `:Max` or `:Min`. Defaults to `:Min`.

 * `cut_oracle::SDDP.AbstractCutOracle`
 The cut oracle is responsible for collecting and storing the cuts that define
 a value function. The cut oracle may decide that only a subset of the total
 discovered cuts are relevant, which improves solution speed by reducing the
 size of the subproblems that need solving. Currently must be one of
    * `DefaultCutOracle()` (see `DefaultCutOracle` for explanation)
    * `LevelOneCutOracle()`(see `LevelOneCutOracle` for explanation)

 * `risk_measure`
 If a single risk measure is given (i.e. `risk_measure = Expectation()`), then
 this measure will be applied to every stage in the problem. Another option is
 to provide a vector of risk measures. There must be one element for every
 stage. For example:

    risk_measure = [ NestedAVaR(lambda=0.5, beta=0.25), Expectation() ]

will apply the `i`'th element of `risk_measure` to every Markov state in the
`i`'th stage. The last option is to provide a vector (one element for each
stage) of vectors of risk measures (one for each Markov state in the stage).
For example:

    risk_measure = [
    # Stage 1 Markov 1 # Stage 1 Markov 2 #
        [ Expectation(), Expectation() ],
        # ------- Stage 2 Markov 1 ------- ## ------- Stage 2 Markov 2 ------- #
        [ NestedAVaR(lambda=0.5, beta=0.25), NestedAVaR(lambda=0.25, beta=0.3) ]
        ]

Note that even though the last stage does not have a future cost function
associated with it (as it has no children), we still have to specify a risk
measure. This is necessary to simplify the implementation of the algorithm.

For more help see `NestedAVaR` or `Expectation`.

 * `markov_transition`
Define the transition probabilties of the stage graph. If a single array is
given, it is assumed that there is an equal number of Markov states in each
stage and the transition probabilities are stage invariant. Row indices
represent the Markov state in the previous stage. Column indices represent the
Markov state in the current stage. Therefore:

    markov_transition = [0.1 0.9; 0.8 0.2]

is the transition matrix when there is 10% chance of transitioning from Markov
state 1 to Markov state 1, a 90% chance of transitioning from Markov state 1
to Markov state 2, an 80% chance of transitioning from Markov state 2 to Markov
state 1, and a 20% chance of transitioning from Markov state 2 to Markov state 2.

# Returns
 * `m`: the `SDDPModel`
"""
function SDDPModel(build!::Function;
    sense                = :Min,
    stages::Int          = 1,
    objective_bound      = nothing,
    markov_transition    = [ones(Float64, (1,1)) for t in 1:stages],
    risk_measure         = Expectation(),
    cut_oracle::AbstractCutOracle = DefaultCutOracle(),
    solver               = UnsetSolver(),
    value_function       = DefaultValueFunction(cut_oracle),
    )
    if objective_bound == nothing
        error("You must specify the objective_bound keyword")
    end
    # check number of arguments to SDDPModel() do [args...] ...  model def ... end
    num_args = n_args(build!)
    if num_args < 2
        error("""Too few arguments in
            SDDPModel() do args...
            end""")
    elseif num_args > 3
        error("""Too many arguments
            SDDPModel() do args...
            end""")
    end

    # New SDDPModel
    m = newSDDPModel(sense, getel(AbstractValueFunction, value_function, 1, 1), build!)

    for t in 1:stages
        markov_transition_matrix = getel(Array{Float64, 2}, markov_transition, t)
        # check that
        if num_args == 2 && size(markov_transition_matrix, 2) > 1
            error("""Because you specified a noise tree in the SDDPModel constructor, you need to use the

                SDDPModel() do sp, stage, markov_state
                    ... model definition ...
                end

            syntax.""")
        end
        stage = Stage(t,markov_transition_matrix)
        for i in 1:size(markov_transition_matrix, 2)
            mod = Subproblem(
                finalstage     = (t == stages),
                stage          = t,
                markov_state   = i,
                sense          = optimisationsense(sense),
                bound          = float(getel(Real, objective_bound, t, i)),
                risk_measure   = getel(AbstractRiskMeasure, risk_measure, t, i),
                value_function = deepcopy(getel(AbstractValueFunction, value_function, t, i))
            )
            setsolver(mod, getel(JuMP.MathProgBase.AbstractMathProgSolver, solver, t, i))
            # dispatch to correct function
            # maybe we should do this with tuples
            if num_args == 3
                build!(mod, t, i)
            else
                build!(mod, t)
            end
            # # Uniform noise probability for now
            if length(ext(mod).noises) != length(ext(mod).noiseprobability)
                setnoiseprobability!(mod, ones(length(ext(mod).noises)) / length(ext(mod).noises))
            end
            push!(stage.subproblems, mod)
        end
        push!(m.stages, stage)
    end
    m
end

function forwardpass!(m::SDDPModel, settings::Settings, solutionstore=nothing)
    last_markov_state = 1
    noiseidx = 0
    obj = 0.0
    for (t, stage) in enumerate(stages(m))
        # choose markov state
        (last_markov_state, sp) = samplesubproblem(stage, last_markov_state, solutionstore)
        if t > 1
            setstates!(m, sp)
        end

        # choose and set RHS noise
        if hasnoises(sp)
            (noiseidx, noise) = samplenoise(sp, solutionstore)
            setnoise!(sp, noise)
        end
        # solve subproblem
        solvesubproblem!(ForwardPass, m, sp, solutionstore)
        # store stage obj
        obj += getstageobjective(sp)
        # store state
        savestates!(stage.state, sp)
        # save solution for simulations (defaults to no-op)
        savesolution!(solutionstore, last_markov_state, noiseidx, sp, t)
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
            @timeit TIMER "Cut addition" begin
                modifyvaluefunction!(m, settings, sp)
            end
        end
    end
    reset!(m.storage)
    for sp in subproblems(m, 1)
        solvesubproblem!(BackwardPass, m, sp)
    end
    return dot(m.storage.objective, m.storage.probability)
end

function iteration!(m::SDDPModel, settings::Settings)
    t = time()
    @timeit TIMER "Forward Pass" begin
        simulation_objective = forwardpass!(m, settings)
    end
    time_forwards = time() - t
    @timeit TIMER "Backward Pass" begin
        objective_bound = backwardpass!(m, settings)
    end
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

function rebuild!(m::SDDPModel)
    for (t, stage) in enumerate(stages(m))
        if t == length(stages(m))
            continue
        end
        for sp in subproblems(stage)
            rebuildsubproblem!(m, sp)
        end
    end
end
function JuMP.solve(::Serial, m::SDDPModel, settings::Settings=Settings())
    status = :solving
    time_simulating, time_cutting = 0.0, 0.0
    objectives = CachedVector(Float64)
    nsimulations, iteration, keep_iterating = 0, 1, true
    start_time = time()
    while keep_iterating
        # add cuts
        @timeit TIMER "Iteration Phase" begin
            (objective_bound, time_backwards, simulation_objective, time_forwards) = iteration!(m, settings)
        end
        # update timers and bounds
        time_cutting += time_backwards + time_forwards
        lower, upper = simulation_objective, simulation_objective

        if applicable(iteration, settings.cut_selection_frequency)
            @timeit TIMER "Cut Selection" begin
                # run cut selection
                settings.print_level > 1 && info("Running Cut Selection")
                rebuild!(m)
            end
        end

        if applicable(iteration, settings.simulation.frequency)
            t = time()
            @timeit TIMER "Simulation Phase" begin
                # simulate policy
                settings.print_level > 1 && info("Running Monte-Carlo Simulation")
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
            end
            time_simulating += time() - t
        end

        total_time = time() - start_time

        addsolutionlog!(m, settings, iteration, objective_bound, lower, upper, time_cutting, nsimulations, time_simulating, total_time, !applicable(iteration, settings.simulation.frequency))

        status, keep_iterating = testboundstall(m, settings, status, keep_iterating)

        if total_time > settings.time_limit
            status = :time_limit
            keep_iterating = false
        end

        iteration += 1
        if iteration > settings.iteration_limit
            status = :iteration_limit
            keep_iterating = false
        end

    end
    status
end

function addsolutionlog!(m, settings, iteration, objective, lower, upper, cutting_time, simulations, simulation_time, total_time, printsingle)
    push!(m.log, SolutionLog(iteration, objective, lower, upper, cutting_time, simulations, simulation_time, total_time))
    print(print, settings, m.log[end], printsingle, m.sense == :Min)
end

function testboundstall(m::SDDPModel, settings::Settings, status::Symbol, keep_iterating::Bool)
    if keep_iterating
        if settings.bound_convergence.iterations > 1 && length(m.log) >= settings.bound_convergence.iterations
            last_n = map(l->l.bound, m.log[end-settings.bound_convergence.iterations+1:end])
            if all(last_n - mean(last_n) .< settings.bound_convergence.atol) || all(abs.(last_n / mean(last_n)-1) .<    settings.bound_convergence.rtol)
                return :bound_convergence, false
            end
        end
    end
    return status, keep_iterating
end

"""
    solve(m::SDDPModel; kwargs...)

# Description

Solve the SDDPModel `m` using SDDP. Accepts a number of keyword arguments to
control the solution process.

# Positional arguments
 * `m`: the SDDPModel to solve

# Keyword arguments
 * `iteration_limit::Int`:
    The maximum number of cuts to add to a single stage problem before terminating.
 * `time_limit::Real`:
    The maximum number of seconds to compute for before termination.
    Defaults to `Inf`.
 * `simulation::MonteCarloSimulation`: see `MonteCarloSimulation`
 * `bound_convergence::BoundConvergence`: see `BoundConvergence`
 * `cut_selection_frequency::Int`:
    Frequency (by iteration) with which to rebuild subproblems using a subset of
    cuts. Frequent cut selection (i.e. `cut_selection_frequency` is small) reduces
    the size of the subproblems that are solved, but incurrs the overhead of rebuilding
    the subproblems. However, infrequent cut selection (i.e.
    `cut_selection_frequency` is large) allows the subproblems to grow large (many
    constraints) leading to an increase in the solve time of individual subproblems.
    Defaults to `0` (never run).
 * `print_level::Int`:
     0 - off: nothing logged.
     1 - on: solve iterations logged.
     2 - verbose: detailed timing information is also logged.
     Defaults to `1`
 * `log_file::String`:
    Relative filename to write the log to disk. Defaults to `""` (no log written)
 * `solve_type`:
    One of
    * `Asynchronous()` - solve using a parallelised algorithm
    * `Serial()` - solve using a serial algorithm
    Default chooses automatically based on the number of available processors.
 * `reduce_memory_footprint::Bool`:
    Implements the idea proposed in https://github.com/JuliaOpt/JuMP.jl/issues/969#issuecomment-282191105
    to reduce the memory consumption when running SDDP. This is an issue if you
    wish to save the model `m` to disk since it discards important information.
    Defaults to `false`.
 * `cut_output_file::String`:
    Relative filename to write discovered cuts to disk. Defaults to `""` (no cuts written)

# Returns
 * `status::Symbol`:
    Reason for termination. One of
    * `:solving`
    * `:interrupted`
    * `:converged`
    * `:iteration_limit`
    * `:bound_convergence`
    * `:time_limit`

"""
function JuMP.solve(m::SDDPModel;
        iteration_limit::Int      = Int(1e9),
        max_iterations::Union{Int, Void} = nothing,
        time_limit::Real          = Inf, # seconds
        simulation = MonteCarloSimulation(
                frequency   = 0,
                min         = 20,
                step        = 1,
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
        print_level::Int             = 1,
        log_file::String             = "",
        # automatically chose Serial or Asynchronous
        solve_type::SDDPSolveType    = nprocs()>2 ? Asynchronous() : Serial(),
        # this reduces memory but you shouldn't use it if you want to save the
        # sddp model since it throws away some information
        reduce_memory_footprint      = false,
        cut_output_file::String      = ""
    )
    if max_iterations != nothing
        warn("The keyword `max_iterations` is deprecated. Use `iteration_limit` instead.")
        iteration_limit = max_iterations
    end
    reset_timer!(TIMER)

    cut_output_file_handle = if cut_output_file != ""
        open(cut_output_file, "w")
    else
        ff = IOStream("")
        close(ff)
        ff
    end

    settings = Settings(
        iteration_limit,
        time_limit,
        simulation,
        bound_convergence,
        cut_selection_frequency,
        print_level,
        log_file,
        reduce_memory_footprint,
        cut_output_file_handle,
        isa(solve_type, Asynchronous)
    )

    print(printheader, settings, m, solve_type)
    status = :solving
    try
        timeit(TIMER, "Solve") do
            status = solve(solve_type, m, settings)
        end
    catch ex
        if isa(ex, InterruptException)
            warn("Terminating solve due to user interaction")
            status = :interrupted
        else
            rethrow(ex)
        end
    finally
        close(cut_output_file_handle)
    end
    print(printfooter, settings, m, settings, status, TIMER)
    status
end

end
