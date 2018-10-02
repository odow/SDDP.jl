function print_banner(io=stdout)
    println(io, "———————————————————————————————————————————————————————————————————————————————")
    println(io, "                         Kokako - © Oscar Dowson, 2018.                        ")
    println(io, " Iteration Simulation Bound")
    println(io, "———————————————————————————————————————————————————————————————————————————————")
end

function print_iteration(iteration, simulation, bound, io=stdout)
    println(io, "$(iteration), $(simulation), $(bound)")
end
