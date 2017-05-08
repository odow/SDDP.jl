#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export Asyncronous

struct Asyncronous <: SDDPSolveType end

function sendtoworkers(;args...)
    for p in workers()
        for (nm, val) in args
            @spawnat(p, eval(SDDP, Expr(:(=), nm, deepcopy(val))))
        end
    end
end

function async_iteration!{C}(T, settings::Settings, slave::Vector{C})
    m = SDDP.m::T
    while length(slave) > 0
        c = pop!(slave)
        addcut!(m, c)
    end
    if !haskey(m.ext, :cuts)
        m.ext[:cuts] = C[]
    else
        empty!(m.ext[:cuts])
    end
    (objective_bound, time_backwards, simulation_objective, time_forwards) = iteration!(m, settings)
    y = similar(m.ext[:cuts])
    copy!(y, m.ext[:cuts])
    empty!(m.ext[:cuts])
    y, objective_bound, simulation_objective
end
function async_forwardpass!(T, settings::Settings)
    m = SDDP.m::T
    forwardpass!(m, settings)
end
function JuMP.solve(::Asyncronous, m::SDDPModel, settings::Settings=Settings())
    status = :solving
    iterationtype = :cutting
    iteration = 1
    simulations = 0
    start_time = time()
    cutting_time = 0.0
    simulation_time = 0.0
    best_objective = Inf
    objectives = CachedVector(Float64)
    simidx = 1

    np = nprocs()
    if np <= 2
        error("You've only loaded two processes. You should solve using the Serial solver instead.")
    end
    slaves = Dict{Int, Vector{Tuple{Int, Int, Cut}}}()
    for i in workers()
        slaves[i] = Tuple{Int, Int, Cut}[]
    end

    begin
        t=time()
        sendtoworkers(m=deepcopy(m))
        info("Took $(round(time() - t, 2)) seconds to copy model to all processes.")
    end

    nextiter() = (nidx=iteration;iteration+=1;nidx)
    getiter() = (iteration)
    # nextsim!() =(simidx+=1)
    addsimulation!() = (simulations += 1; s=simulations;s)
    settype!(x) = (iterationtype = x)
    setbestobjective!(v) = (best_objective = v)
    cutting_timer = -1.0
    simulation_timer = -1.0
    @sync begin
        for p in workers()
            @async begin
                while true
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
                        (cuts, objective_bound, simulation_objective) = remotecall_fetch(async_iteration!, p, typeof(m), settings, newcuts)
                        for cut in cuts
                            addcut!(m, cut)
                            for p2 in workers()
                                p == p2 && continue
                                push!(slaves[p2], cut)
                            end
                        end

                        total_time = time() - start_time
                        if objective_bound < best_objective
                            setbestobjective!(objective_bound)
                        end
                        cutting_time += time() - cutting_timer
                        addsolutionlog!(m, settings, it, best_objective, simulation_objective, simulation_objective, cutting_time , simulations, simulation_time, total_time, true)

                        status, keep_iterating = testconvergence(m, settings)
                        !keep_iterating && break


                    elseif iterationtype == :simulation

                        if simidx > length(settings.simulation.steps)
                            settype!(:cutting)
                            continue
                        end
                        simulation_timer = time()

                        push!(objectives, remotecall_fetch(async_forwardpass!, p, typeof(m), settings))
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

# addprocs(4)
#
# @everywhere begin
#
#     struct Cut
#         x::Float64
#         p::Int
#     end
#
#     struct Model
#         x::Vector{Cut}
#     end
#     Model() = Model(Cut[])
#
#     function foo(m::Model, slave::Channel{Cut})
#         while isready(slave)
#             push!(m.x, take!(slave))
#         end
#         cut = Cut(rand(), myid())
#         sleep(1)
#         push!(m.x, cut)
#         cut
#     end
#     function mainfoo(slave)
#         m = Main.m::Model
#         foo(m, slave)
#     end
# end
#
#
# function mastermanager!(m::Model, N)
#     np = nprocs()
#     slaves = [Channel{Cut}(2 * np) for p in 1:np]
#     sendto(workers(), m=deepcopy(m))
#     n = 1
#     nextcut() = (nidx=n;n+=1;nidx)
#     @sync begin
#         for p in 1:np
#             p == myid() && continue
#             @async begin
#                 for i in 1:4
#                     idx = nextcut()
#                     if idx > N
#                         break
#                     end
#                     cut = remotecall_fetch(mainfoo, p, slaves[p])
#                     println("discovered $cut")
#                     push!(m.x, cut)
#                     for p2 in 1:np
#                         cut.p == p2 && continue
#                         put!(slaves[p2], cut)
#                     end
#                 end
#             end
#         end
#     end
#     for s in slaves
#         close(s)
#     end
# end
#
# m = Model()
# @time mastermanager!(m, 10)
# @assert length(m.x) == 10
