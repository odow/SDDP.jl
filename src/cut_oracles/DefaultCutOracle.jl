#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

"""
    DefaultCutOracle()

# Description

Initialize the default cut oracle.

This oracle keeps every cut discovered and does not perform cut selection.
"""
struct DefaultCutOracle <: AbstractCutOracle
    cuts::Vector{Cut}
end
DefaultCutOracle() = DefaultCutOracle(Cut[])

function store_cut(oracle::DefaultCutOracle, model::SDDPModel,
                   subproblem::JuMP.Model, cut::Cut)
    push!(oracle.cuts, cut)
    return
end
valid_cuts(oracle::DefaultCutOracle) = oracle.cuts
all_cuts(oracle::DefaultCutOracle) = oracle.cuts
