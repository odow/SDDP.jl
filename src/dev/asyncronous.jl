#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

master node just stores a jagged array of cuts
    cuts[stage][markovstate][pricestate] = Cut[]

slave nodes run independent SDDPModels.
    getcuts!(m::SDDPModel, stage, markovstate, pricestate, n::Int) # n is current number
    passcuts!(m::SDDPModel, stage, markovstate, pricestate)


asyncronous solve goes

while not converged

    @async do a whole lot of cuts

    @async simulate bound

    check convergence
end
