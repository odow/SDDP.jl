#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

function printheader(io::IO)
    println(io, "               Objective               |  Cut  Passes  |  Simulations  ")
    println(io, "      Expected      |   Bound   % Gap  |   #     Time  |    #    Time  ")
end

function Base.print(io::IO, l::SolutionLog, printmean::Bool=false)
    if printmean
        bound_string = string("     ", humanize(0.5 * (l.lower_statistical_bound + l.upper_statistical_bound), "8.3f"), "     ")
        rtol_string = "      "
    else
        bound_string = string(
            humanize(l.lower_statistical_bound, "8.3f"), " ",
            humanize(l.upper_statistical_bound, "8.3f")
        )
        rtol_string = humanize(100*rtol(l.bound, l.upper_statistical_bound), "5.1f")
    end

    println(io,
        @sprintf("%s | %s %s | %s %s | %s %s",
            bound_string,
            humanize(l.bound, "8.3f"),
            rtol_string,
            humanize(l.iteration),
            humanize(l.timecuts),
            humanize(l.simulations),
            humanize(l.timesimulations)
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

function printheader(s::String)
    open(s, "a") do file
        printheader(file)
    end
end
function Base.print(s::String, l::SolutionLog, printmean::Bool=false)
    open(s, "a") do file
        print(file, l, printmean)
    end
end
