#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

struct DefaultValueFunction{C<:AbstractCutOracle} <: AbstractValueFunction
    cutmanager::C
    theta::JuMP.Variable
end

function DefaultValueFunction(m::JuMP.Model, sense, bound, cutmanager=DefaultCutOracle())
    DefaultValueFunction(cutmanager, futureobjective!(sense, m, bound))
end

function init!(vf::Type{DefaultValueFunction}, m::JuMP.Model, sense, bound, cutmanager)
    DefaultValueFunction(m::JuMP.Model, sense, bound, cutmanager)
end
