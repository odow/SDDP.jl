#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

stageobjective!(vf::AbstractValueFunction, sp::JuMP.Model, obj...) = error("You need this method")
getstageobjective(vf::AbstractValueFunction, sp::JuMP.Model) = error("You need this method")
init!(vf::AbstractValueFunction, m::JuMP.Model, sense, bound, cutmanager) = error("You need this method")
modifyvaluefunction!(vf::AbstractValueFunction, m::SDDPModel, sp::JuMP.Model) = error("You need this method")
rebuildsubproblem!(vf::AbstractValueFunction, m::SDDPModel, sp::JuMP.Model) = nothing

stageobjective!(sp::JuMP.Model, obj...) = stageobjective!(valueoracle(sp), sp, obj...)
getstageobjective(sp::JuMP.Model) = getstageobjective(valueoracle(sp), sp)
modifyvaluefunction!(m::SDDPModel, sp::JuMP.Model) = modifyvaluefunction!(valueoracle(sp), m, sp)
rebuildsubproblem!(m::SDDPModel, sp::JuMP.Model) = rebuildsubproblem!(valueoracle(sp), m, sp)

function constructcut(m::SDDPModel, sp::JuMP.Model)
    # theta <=/>= E[ (y - πᵀx̄) + πᵀx ]
    intercept = 0.0
    coefficients = zeros(nstates(sp))
    for i in 1:length(m.storage.objective)
        intercept += m.storage.modifiedprobability[i] * (m.storage.objective[i] - dot(getstage(m, ext(sp).stage).state, m.storage.duals[i]))
        for j in 1:nstates(sp)
            coefficients[j] += m.storage.modifiedprobability[i] * m.storage.duals[i][j]
        end
    end
    Cut(intercept, coefficients)
end

_addcut!(::Type{Min}, sp, theta, affexpr) = @constraint(sp, theta >= affexpr)
_addcut!(::Type{Max}, sp, theta, affexpr) = @constraint(sp, theta <= affexpr)
