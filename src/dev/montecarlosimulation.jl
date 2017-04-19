#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

function simulate(m::SDDPModel, n::Int, variablestorecord::Vector{Symbol}=Symbol[])
    for stage 1 to t-1
        solve stage i
        pass new state forward
        pass new price forward
        transition to new markov state
        record variables
        record stage cost
    end

    solve stage t
    record variables
    record stage cost
end
