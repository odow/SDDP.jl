#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

"""
    store_cut(oracle::AbstactCutOracle, model::SDDPModel,
              subproblem::JuMP.Model, cut::Cut)

# Description

Store the cut `cut` in the Cut Oracle `oracle`. `oracle` will belong to
`subproblem` in the SDDPModel `model`.
"""
function store_cut end

"""
    valid_cuts(oracle::AbstactCutOracle)

# Description

Return an iterable list of all the valid cuts contained within `oracle`.
"""
function valid_cuts end

"""
    all_cuts(oracle::AbstactCutOracle)

# Description

Return an iterable list of *all* the cuts contained within `oracle`, not just
those that are returned by `valid_cuts`.
"""
function all_cuts end

include("DefaultCutOracle.jl")
include("LevelOneCutOracle.jl")

@deprecate storecut! store_cut
@deprecate validcuts valid_cuts
@deprecate allcuts all_cuts
