#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

"""
    storecut!(oracle::AbstractCutOracle, cut::Cut)

    This function adds the cut to the CutOracle.
"""
storecut!(oracle::AbstractCutOracle, cut::Cut) = error("""
    You must define the function
        storecut!(oracle::$(typeof(oracle)), m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int, cut)
    that is overloaded for your oracle of type $(typeof(oracle)).
""")
# more expansive method that can also be overloaded
storecut!(oracle::AbstractCutOracle, m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int, cut) = storecut!(oracle, cut)

"""
    validcuts(oracle::AbstactCutOracle)
    This function returns an iterable list of all the valid cuts contained within the oracle.
"""
validcuts(oracle::AbstractCutOracle) = error("""You must define the function validcuts(oracle::$(typeof(oracle))) that is overloaded for your
    oracle of type $(typeof(oracle)).""")


struct DefaultCutOracle <: AbstractCutOracle
    cuts::Vector{Cut}
end
DefaultCutOracle() = DefaultCutOracle(Cut[])

storecut!(oracle::DefaultCutOracle, cut::Cut) = push!(oracle.cuts, cut)
validcuts(oracle::DefaultCutOracle) = oracle.cuts
