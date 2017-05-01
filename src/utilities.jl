#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
getstage(m::SDDPModel, t) = stages(m)[t]

subproblems(s::Stage) = s.subproblems
getsubproblem(s::Stage, i::Int) = subproblems(s)[i]
subproblems(m::SDDPModel, t) = subproblems(getstage(m, t))
getsubproblem(m::SDDPModel, t, i) = getsubproblem(getstage(m, t), i)

nstages(m::SDDPModel) = length(stages(m))
nsubproblems(m::SDDPModel, t::Int) = length(subproblems(m, t))

function n_args(f::Function)
    @assert length(methods(f)) == 1
    return methods(f).mt.max_args-1
end
getel{A, T <: A}(::Type{A}, x::T, t::Int, i::Int) = x
getel{A, T <: A}(::Type{A}, x::Vector{T}, t::Int, i::Int) = x[t]
getel{A, T <: A}(::Type{A}, x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]
getel{A, T <: A}(::Type{A}, x::T, t::Int) = x
getel{A, T <: A}(::Type{A}, x::Vector{T}, t::Int) = x[t]

# Construct a confidence interval
function confidenceinterval(x, conf_level=0.95)
    tstar = quantile(TDist(length(x)-1), 1 - (1 - conf_level)/2)
    err = tstar * std(x)/sqrt(length(x))
    mu = mean(x)
    return mu - err, mu + err
end

function rtol(x, y)
    if abs(y) < 1e-6
        return x - y
    else
        (x - y) / y
    end
end

function sample(x::AbstractVector{Float64})
    r = rand()
    for i in 1:length(x)
        @inbounds r -= x[i]
        if r < eps(Float64)
            return i
        end
    end
    error("x must be a discrete probablity distribution that sums to one. sum= $(sum(x))")
end
