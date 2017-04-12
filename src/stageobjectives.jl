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

addfutureobjective!(::Type{Max}, m::JuMP.Model) = @variable(m, upperbound = ext(m).problembound)
addfutureobjective!(::Type{Min}, m::JuMP.Model) = @variable(m, lowerbound = ext(m).problembound)
addfutureobjective!(m::JuMP.Model) = addfutureobjective!(ext(m).sense, m)

function stageobjective!(sp::JuMP.Model, obj::AffExpr)
    @assert length(ext(sp).valuefunctions) == 0
    push!(ext(sp).valuefunctions, ValueFunction(0.0, addfutureobjective!(sp)))
    JuMP.setobjective(sp, getsense(sp), obj)
end
stageobjective!(sp::JuMP.Model, obj::JuMP.Variable) = stageobjective!(sp, AffExpr(obj))

function stageobjective!(sp::JuMP.Model,        # the JuMP model
                    ribs::Vector{Float64},      # locations of ribs
                    pricetransition::Function,  # (old price, noise) -> (new price)
                    noises::AbstractVector,     # vector of noises
                    objectivefunction::Function # R -> AffExpr
                )
    @assert length(ext(sp).valuefunctions) == 0
    sort!(ribs)
    for r in ribs
        push!(ext(sp).valuefunctions, ValueFunction(r, addfutureobjective!(sp)))
    end
    priceoracle = ext(sp).priceoracle
    priceoracle.pricetransition = transition
    priceoracle.pricescenarios = noises
    priceoracle.objective = objectivefunction
end
