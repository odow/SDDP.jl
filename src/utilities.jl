#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################
getsense(::Type{Max}) = :Max
getsense(::Type{Min}) = :Min
getsense(m::JuMP.Model) = getsense(ext(m).sense)
function optimisationsense(s::Symbol)
    if s==:Min
        return Min
    elseif s==:Max
        return Max
    else
        error("Unknown optimisation sense $s. Must be :Max or :Min")
    end
end

futureobjective!(::Type{Max}, m::JuMP.Model, bound) = @variable(m, upperbound = bound)
futureobjective!(::Type{Min}, m::JuMP.Model, bound) = @variable(m, lowerbound = bound)

stages(m::SDDPModel) = m.stages
stage(m::SDDPModel, t) = stages(m)[t]

subproblems(m::SDDPModel, t) = stage(m, t).subproblems
subproblem(m::SDDPModel, t, i) = subproblem(m, t)[i]

nstages(m::SDDPModel) = length(stages(m))
nsubproblems(m::SDDPModel, t::Int) = length(subproblems(m, t))
