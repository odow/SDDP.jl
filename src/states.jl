#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

const comparison_symbols = [:(<=), :(>=), :(==)]
is_comparison(x) = Base.Meta.isexpr(x, :comparison) || (Base.Meta.isexpr(x, :call) && x.args[1] in comparison_symbols)

setvalue!(x::State, y::Float64) = JuMP.setRHS(x.constraint, y)
JuMP.getdual(x::State) = JuMP.getdual(x.constraint)
states(sp::JuMP.Model) = sp.ext[:SDDP].states

function statevariable!(m::JuMP.Model, xin::JuMP.Variable, xout::JuMP.Variable)
    push!(states(m),
            State(
                xout,
                @constraint(m, xin == getvalue(xout))
            )
    )
end

function statevariable!(m::JuMP.Model, xin::Array{JuMP.Variable}, xout::Array{JuMP.Variable})
    @assert length(xin) == length(xout)
    for i in 1:length(xin)
        statevariable!(m, xin[i], xout[i])
    end
end

function statevariable!{T<:Union{JuMP.JuMPArray, JuMP.JuMPDict}}(sp::JuMP.Model, xin::T, xout::T)
    @assert length(keys(xin)) == length(keys(xout))
    for key in keys(xin)
        statevariable!(sp, xin[key...], xout[key...])
    end
end

"""
    @state(sp, stateleaving, stateentering)
Define a new state variable in the subproblem `sp`.
Arguments:
    sp               the subproblem
    stateleaving     any valid JuMP `@variable` syntax to define the value of the state variable at the end of the stage
    stateentering    any valid JuMP `@variable` syntax to define the value of the state variable at the beginning of the stage
Usage:
    @state(sp, 0 <= x[i=1:3] <= 1, x0=rand(3)[i] )
    @state(sp,      y        <= 1, y0=0.5        )
    @state(sp,      z            , z0=0.5        )
"""
macro state(sp, x, x0)
    sp = esc(sp)                        # escape the model
    @assert x0.head == :(=)             # must be a keyword
    symin, rhs = x0.args                # name of the statein variable
    if is_comparison(x)
        if length(x.args) == 5          # double sided
            xin = identity(x.args[3])       # variable is in middle
        elseif length(x.args) == 3      # single comparison
            xin = identity(x.args[2])       # variable is second entry
        else
            error("Unknown format for $(x)")
        end
    else
        xin = identity(x)                   # no bounds
    end
    if isa(xin, Expr)                   # x has indices
        xin.args[1] = symin             # so just change the name
    else                                # its just a Symbol
        xin = symin                     # so change the Symbol
    end
    # ex = Expr(:(=), QuoteNode(:start), esc(rhs))
    # @show ex, dump(ex)
    quote
        stateout = $(Expr(:macrocall, Symbol("@variable"), sp, esc(x), :(start=esc(rhs))))
        statein  = $(Expr(:macrocall, Symbol("@variable"), sp, esc(xin)))
        statevariable!($sp, statein, stateout)
        stateout, statein
    end
end
