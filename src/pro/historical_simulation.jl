#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export historicalsimulation

function historicalsimulation{C}(m::SDDPModel{DefaultValueFunction{C}},
        variables::Vector{Symbol} = Symbol[];
        scenarios::Vector{Int}    = ones(Int, length(m.stages)),
        markovstates::Vector{Int} = ones(Int, length(m.stages))
    )

    store = newsolutionstore(variables)

    obj = 0.0

    for t in 1:length(m.stages)
        sp = getsubproblem(m, t, markovstates[t])
        setscenario!(sp, ext(sp).scenarios[scenarios[t]])
        solvesubproblem!(ForwardPass, m, sp)
        # solve subproblem
        solvesubproblem!(ForwardPass, m, sp)
        # store stage obj
        obj += getstageobjective(sp)
        # store state
        savestates!(getstage(m, t).state, sp)
        # save solution for simulations (defaults to no-op)
        savesolution!(store, markovstates[t], scenarios[t], sp)
    end
    store[:objective] = obj
    return store
end
