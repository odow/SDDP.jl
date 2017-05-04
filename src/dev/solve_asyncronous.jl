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
            @spawnat(p, eval(Main, Expr(:(=), nm, val)))
        end
    end
end

function async_iteration!{C}(T, settings::Settings, slave::Channel{C})
    m = SDDP.m::T
    while isready(slave)
        c = take!(slave)
        addcut!(m, c)
    end
    (objective_bound, time_backwards, simulation_objective, time_forwards) = iteration!(m, settings)

end

function solve(::Asyncronous, m::SDDPModel, settings::Settings=Settings())
    status = :solving
    time_simulating, time_cutting = 0.0, 0.0
    objectives = CachedVector(Float64)
    nsimulations, iteration, keep_iterating = 0, 1, true
    start_time = time()

    np = nprocs()

    if np <= 2
        error("You've only loaded two processes. You should solve using the Serial solver instead.")
    end

    nextiter() = (nidx=iteration;iteration+=1;nidx)
    @sync begin
        for p in workers()
            @async begin
                while true
                    it = nextiter()
                    if it > settings.maximum_iterations
                        break
                    end
                    asyncreturn = remotecall_fetch(async_iteration!, p, typeof(m), settings, slaves[p])
                    # push!(m.x, cut)
                    for p2 in workers()
                        cut.p == p2 && continue
                        # put!(slaves[p2], cut)
                    end
                end
            end
        end
    end

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



println("The end")
