#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

function printheader(io::IO)
    println(io, "               Objective               |  Cut  Passes  |  Simulations  ")
    println(io, "      Expected      |   Bound   % Gap  |  Cuts   Time  |  Sims   Time  ")
end

function Base.print(io::IO, l::SolutionLog)
    if l.lower_statistical_bound == l.upper_statistical_bound
        bound_string = string("     ", humanize(l.lower_statistical_bound, "8.3f"), "     ")
    else
        bound_string = string(
            humanize(l.lower_statistical_bound, "8.3f"), " ",
            humanize(l.upper_statistical_bound, "8.3f")
        )
    end

    println(io,
        @sprintf("%s | %s %s | %s %s | %s %s",
            bound_string,
            humanize(l.bound, "8.3f"),
            "      ",
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

# printheader(STDOUT)
# print(STDOUT, SolutionLog())


# getCloseCIBound(::Type{Min}, m::SDDPModel) = m.confidence_interval[1]
# getCloseCIBound(::Type{Max}, m::SDDPModel) = m.confidence_interval[2]
# getCloseCIBound{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}) = getCloseCIBound(X, m)
# getBound(m::SDDPModel) = m.valid_bound
# atol(::Type{Min}, m::SDDPModel) = getCloseCIBound(m) - getBound(m)
# atol(::Type{Max}, m::SDDPModel) = getBound(m) - getCloseCIBound(m)
# atol{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}) = atol(X, m)
#
# """
#     rtol(SDDPModel)
# Relative tolerance of the solution
# Defined as [Outer bound - closest simulated bound] / [Outer bound]
# """
function rtol(l::SolutionLog)
    abs(getBound(m)) != Inf?atol(m) / abs(getBound(m)):Inf
end
