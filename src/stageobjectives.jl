#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

function stageobjective!{S,V<:DefaultValueFunction,R}(ex::SubproblemExt{S,V,R}, sp::JuMP.Model, obj::AffExpr)
    JuMP.setobjective(sp, getsense(sp), obj + ex.valueoracle.theta)
end

stageobjective!(sp::JuMP.Model, obj::AffExpr) = stageobjective!(ext(sp), sp, obj)
stageobjective!(sp::JuMP.Model, obj::JuMP.Variable) = stageobjective!(sp, AffExpr(obj))

# function stageobjective!(sp::JuMP.Model,        # the JuMP model
#                     ribs::Vector{Float64},      # locations of ribs
#                     pricetransition::Function,  # (old price, noise) -> (new price)
#                     noises::AbstractVector,     # vector of noises
#                     objectivefunction::Function # R -> AffExpr
#                 )
#     @assert length(ext(sp).valuefunctions) == 0
#     sort!(ribs)
#     for r in ribs
#         push!(ext(sp).valuefunctions, ValueFunction(r, addfutureobjective!(sp)))
#     end
#     priceoracle = ext(sp).priceoracle
#     priceoracle.pricetransition = transition
#     priceoracle.pricescenarios = noises
#     priceoracle.objective = objectivefunction
# end
