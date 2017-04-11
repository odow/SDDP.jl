#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

for stage 1 to t
    solve stage i
    pass new state forward
    pass new price forward
    transition to new markov state
    record state
end

for stage t-1 to stage 1
    for each price rib above and below last price point
        for each new markov state
            for each new price outcome
                for each new scenario outcome
                    solve problem
                    record objective
                    record duals on incoming state variables
                end
            end
        end
        add cut
    end
end
