#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
"""
    Asyncrounous(; kwargs...)

# Define

Type used to dispatch and control the behaviour of the asyncronous solution
algorithm.

# Arguments

 * `slaves::Vector{Int}` the pid's of the slave processes. Defaults to `workers()`
 * `step::Float64` the number of iterations to complete before adding another
    slave. Used to replicated the scenario incrementation behaviour of
    V. de Matos,A. Philpott, E. Finardi, Improving the performance of Stochastic
    Dual Dynamic Programming, Journal of Computational and Applied Mathematics
    290 (2015) 196–208.

# Examples

    Asyncronous() # load on all workers
    Asyncronous(slaves=[2,3,4]) # load slaves on processes 2, 3, and 4
    Asyncronous(step=10) # perform 10 iterations before adding new slave

"""
immutable Asyncronous <: SDDPSolveType
    slaves::Vector{Int} # pid of slave processors
    step::Float64       # number of iterations before introducing another slave
end
function Asyncronous(;slaves=workers(), step=1/(1 + length(slaves)))
    sl = Int[]
    for s in slaves
        if s == myid()
            warn("Current process passed to Asyncronous() as slave. Ignoring.")
        else
            push!(sl, s)
        end
    end
    Asyncronous(sl, step)
end
Base.show(io::IO, async::Asyncronous) = print(io, "Asyncronous solver with $(length(async.slaves)) slave processors and a step of $(async.step)")

function sendto(procs;args...)
    for p in procs
        p == myid() && continue
        for (nm, val) in args
            io = IOBuffer()
            serialize(io, val)
            @spawnat(p, begin
                eval(SDDP, Expr(:(=), :io, io))
                seekstart(SDDP.io)
                eval(SDDP, Expr(:(=), nm, deserialize(SDDP.io)))
                close(SDDP.io)
            end)
            close(io)
        end
    end
end

function async_iteration!{C}(T, settings::Settings, slave::Vector{C})
    m = SDDP.m::T
    async_iteration!(m, settings, slave)
end
function async_iteration!{C}(m::SDDPModel, settings::Settings, slave::Vector{C})
    while length(slave) > 0
        c = pop!(slave)
        addcut!(m, c)
    end
    if !haskey(m.ext, :cuts)
        m.ext[:cuts] = C[]
    else
        empty!(m.ext[:cuts])
    end
    (objective_bound, time_backwards, simulation_objective, time_forwards) = iteration!(m, settings, false)
    y = copy(m.ext[:cuts])
    empty!(m.ext[:cuts])
    y, objective_bound, simulation_objective
end
function async_forwardpass!(T, settings::Settings)
    m = SDDP.m::T
    forwardpass!(m, settings)
end

function rebuild!(T)
    mm = SDDP.m::T
    rebuild!(mm)
end

function JuMP.solve{T}(async::Asyncronous, m::SDDPModel{T}, settings::Settings=Settings())
    status = :solving
    iterationtype = :cutting
    iteration = 1
    simulations = 0
    start_time = time()
    cutting_time = 0.0
    simulation_time = 0.0
    best_objective = -worstcase(m.sense)
    objectives = CachedVector(Float64)
    simidx = 1

    np = length(async.slaves) + 1
    if np <= 2
        error("You've only loaded one slave process. You should solve using the Serial solver instead.")
    end
    storage_type = getcutstoragetype(T)
    slaves = Dict{Int, Vector{storage_type}}()
    for i in workers()
        slaves[i] = storage_type[]
    end

    begin
        @timeit TIMER "Remote Initialization" begin
            sendto(async.slaves, m=m)
        end
    end


    needs_rebuilding = fill(false, np)
    iterations_since_cut_selection = 0

    nextiter() = (nidx=iteration;iteration+=1;nidx)
    getiter() = (iteration)
    # nextsim!() =(simidx+=1)
    addsimulation!() = (simulations += 1; s=simulations;s)
    settype!(x) = (iterationtype = x)
    setbestobjective!(v) = (best_objective = v)
    cutting_timer = -1.0
    simulation_timer = -1.0
    @sync begin
        for p in async.slaves
            p == myid() && continue
            @async begin
                while true
                    # test time limit
                    if time() - start_time > settings.time_limit
                        status = :time_limit
                        break
                    end
                    # check if slave should be used
                    if !(p in async.slaves[1:min(end, ceil(Int, getiter() / async.step))])
                        sleep(1.0)
                        continue
                    end

                    if iterationtype == :cutting
                        cutting_timer = time()
                        it = nextiter()
                        if applicable(it, settings.simulation.frequency)
                            # the next worker should start simulating
                            settype!(:simulation)
                            reset!(objectives)
                            simidx = 1
                        end
                        if it > settings.max_iterations
                            status = :max_iterations
                            break
                        end
                        newcuts = deepcopy(slaves[p])
                        empty!(slaves[p])
                        @timeit TIMER "Iteration Phase" begin
                            (cuts, objective_bound, simulation_objective) = remotecall_fetch(async_iteration!, p, typeof(m), settings, newcuts)
                        end
                        for cut in cuts
                            if isopen(settings.cut_output_file)
                                writeaynccut!(settings.cut_output_file, cut)
                            end
                            addcut!(m, cut)
                            for p2 in async.slaves
                                p == p2 && continue
                                push!(slaves[p2], cut)
                            end
                        end

                        total_time = time() - start_time
                        if dominates(m.sense, objective_bound, best_objective)
                            setbestobjective!(objective_bound)
                        end
                        cutting_time += time() - cutting_timer

                        iterations_since_cut_selection += 1
                        if applicable(iterations_since_cut_selection, settings.cut_selection_frequency)
                            needs_rebuilding .= true
                            iterations_since_cut_selection = 0
                        end
                        if needs_rebuilding[p]
                            @timeit TIMER "Cut Selection" begin
                                remotecall_fetch(rebuild!, p, typeof(m))
                            end
                            needs_rebuilding[p] = false
                        end

                        addsolutionlog!(m, settings, it, best_objective, simulation_objective, simulation_objective, cutting_time , simulations, simulation_time, total_time, true)

                        status, keep_iterating = testboundstall(m, settings, status, true)
                        !keep_iterating && break

                    elseif iterationtype == :simulation

                        if simidx > length(settings.simulation.steps)
                            settype!(:cutting)
                            continue
                        end
                        simulation_timer = time()

                        @timeit TIMER "Simulation Phase" begin
                            push!(objectives, remotecall_fetch(async_forwardpass!, p, typeof(m), settings))
                        end
                        nsimulations = addsimulation!()

                        simulation_time += time() - simulation_timer
                        if length(objectives) >= settings.simulation.steps[min(end,simidx)]
                            simidx += 1
                            (lower, upper) = confidenceinterval(objectives, settings.simulation.confidence)
                            total_time = time() - start_time
                            if contains(best_objective, lower, upper)
                                if length(objectives) >= settings.simulation.steps[end]
                                    addsolutionlog!(m, settings, iteration-1, best_objective, lower, upper, cutting_time, nsimulations, simulation_time, total_time, false)
                                    # terminate with statistical bound convergence
                                    if settings.simulation.termination
                                        status = :converged
                                        break
                                    end
                                end
                            else
                                addsolutionlog!(m, settings, iteration-1, best_objective, lower, upper, cutting_time, nsimulations, simulation_time, total_time, false)
                                settype!(:cutting)
                            end
                        end
                    end

                end
            end
        end
    end
    status
end


getcutstoragetype{C}(::Type{DefaultValueFunction{C}}) = Tuple{Int, Int, Cut}
function addcut!{C}(m::SDDPModel{DefaultValueFunction{C}}, cut::Tuple{Int, Int, Cut})
    sp = getsubproblem(m, cut[1], cut[2])
    vf = valueoracle(sp)
    storecut!(vf.cutmanager, m, sp, cut[3])
    addcut!(vf, sp, cut[3])
end
function storecut!{C}(m::SDDPModel{DefaultValueFunction{C}}, sp::JuMP.Model, cut::Cut)
    if !haskey(m.ext, :cuts)
        m.ext[:cuts] = Tuple{Int, Int, Cut}[]
    end
    push!(m.ext[:cuts], (ext(sp).stage, ext(sp).markovstate, cut))
end

function writeaynccut!(file, cut::Tuple{Int, Int, Cut})
    writecut!(file, cut[3], cut[1], cut[2])
end
