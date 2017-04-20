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
    stageobjective::QuadExpr
    theta::JuMP.Variable
end

function DefaultValueFunction(m::JuMP.Model, sense, bound, cutmanager=DefaultCutOracle())
    DefaultValueFunction(cutmanager, QuadExpr(0), futureobjective!(sense, m, bound))
end

function init!(::Type{DefaultValueFunction}, m::JuMP.Model, sense, bound, cutmanager)
    DefaultValueFunction(m::JuMP.Model, sense, bound, cutmanager)
end

function stageobjective!(::Type{DefaultValueFunction}, sp::JuMP.Model, obj::AffExpr)
    ext(sp).valueoracle.stageobjective += obj
    JuMP.setobjective(sp, getsense(sp), obj + ex.valueoracle.theta)

getstageobjective(::Type{DefaultValueFunction}, sp::JuMP.Model) = getvalue(ext(sp).valueoracle.stageobjective)

stageobjective!{T<:AbstractValueFunction}(::Type{T}, sp::JuMP.Model, obj::AffExpr) = error("You need this method")
getstageobjective{T<:AbstractValueFunction}(::Type{T}, , sp::JuMP.Model) = error("You need this method")
init!{T<:AbstractValueFunction}(::Type{T}, m::JuMP.Model, sense, bound, cutmanager) = error("You need this method")

stageobjective!(sp::JuMP.Model, obj::AffExpr) = stageobjective!(vftype(sp), sp, obj)
stageobjective!(sp::JuMP.Model, obj::JuMP.Variable) = stageobjective!(sp, AffExpr(obj))
getstageobjective(sp::JuMP.Model) = getstageobjective(ext(sp), sp)
