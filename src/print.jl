#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

function printheader{T}(io::IO, m::SDDPModel{T}, solve_type)
    n = length(m.stages)
    println(io, """--------------------------------------------------------------------------------
                          SDDP Solver. © Oscar Dowson, 2017.
    --------------------------------------------------------------------------------""")
    println(io, """    Solver:
            $(solve_type)
        Model:
            Stages:         $(length(m.stages))
            States:         $(nstates(getsubproblem(m, 1, 1)))
            Subproblems:    $(sum(length(subproblems(s)) for s in stages(m)))
            Value Function: $(summarise(T))
    --------------------------------------------------------------------------------""")
    println(io, "              Objective              |  Cut  Passes    Simulations   Total    ")
    println(io, "     Simulation       Bound   % Gap  |   #     Time     #    Time    Time     ")
    println(io, "--------------------------------------------------------------------------------")
end

function printfooter(io::IO, m::SDDPModel, settings, status, timer)
    println(io, "--------------------------------------------------------------------------------")
    if settings.print_level > 1
        print_timer(io, timer, title="Timing statistics")
        print(io, "\n")
    end
    println(io, """    Other Statistics:
            Iterations:         $(m.log[end].iteration)
            Termination Status: $(status)
    ================================================================================""")
end

function Base.print(io::IO, l::SolutionLog, printmean::Bool=false, is_min=true)
    if printmean
        bound_string = string("     ", humanize(0.5 * (l.lower_statistical_bound + l.upper_statistical_bound), "8.3f"), "     ")
        rtol_string = "      "
    else
        bound_string = string(
            humanize(l.lower_statistical_bound, "8.3f"), " ",
            humanize(l.upper_statistical_bound, "8.3f")
        )
        if is_min
            tol = 100*rtol(l.lower_statistical_bound, l.bound)
        else
            tol = -100*rtol(l.upper_statistical_bound, l.bound)
        end
        rtol_string = humanize(tol, "5.1f")
    end

    println(io,
        @sprintf("%s %s %s | %s %s %s %s %s",
            bound_string,
            humanize(l.bound, "8.3f"),
            rtol_string,
            humanize(l.iteration),
            humanize(l.timecuts),
            humanize(l.simulations),
            humanize(l.timesimulations),
            humanize(l.timetotal)
        )
    )
end

function humanize(value::Int)
    if value < 1000 && value > -1000
        return humanize(value, "5d")
    else
        return humanize(value, "5.1f")
    end
end

function printtofile(foo::Function, filename::String, args...)
    open(filename, "a") do file
        foo(file, args...)
    end
end

function Base.print(printfunc::Function, settings::Settings, args...)
    settings.print_level > 0 && printfunc(STDOUT, args...)
    settings.log_file != "" && printtofile(printfunc, settings.log_file, args...)
end
