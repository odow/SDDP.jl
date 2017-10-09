#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

const comparison_symbols = [:(<=), :(>=), :(==)]
is_comparison(x) = Base.Meta.isexpr(x, :comparison) || (Base.Meta.isexpr(x, :call) && x.args[1] in comparison_symbols)

setvalue!(x::State, y::Float64) = JuMP.setRHS(x.constraint, y)
JuMP.getvalue(x::State) = JuMP.getvalue(x.variable)
JuMP.getdual(x::State) = JuMP.getdual(x.constraint)
states(sp::JuMP.Model) = ext(sp).states
nstates(sp::JuMP.Model) = length(states(sp))
function saveduals!(y, sp::JuMP.Model)
    for (i, state) in enumerate(states(sp))
        y[i] = getdual(state)
    end
end
function savestates!(y, sp::JuMP.Model)
    padvec!(y, nstates(sp))
    for (i, state) in enumerate(states(sp))
        y[i] = getvalue(state)
    end
end

function setstates!(m, sp)
    s = getstage(m, ext(sp).stage-1)
    for (st, v) in zip(states(sp), s.state)
        setvalue!(st, v)
    end
end

function padvec!{T}(x::AbstractVector{T}, n::Int)
    append!(x, zeros(T, max(0, n - length(x))))
end

function statevariable!(m::JuMP.Model, xin::JuMP.Variable, xout::JuMP.Variable)
    push!(states(m),
            State(
                xout,
                @constraint(m, xin == getvalue(xin))
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

_copy(x::Symbol) = x
_copy(x::Expr) = copy(x)

"""
    @state(sp, stateleaving, stateentering)

# Description

Define a new state variable in the subproblem `sp`.

# Arguments
 * `sp`               the subproblem
 * `stateleaving`     any valid JuMP `@variable` syntax to define the value of the state variable at the end of the stage
 * `stateentering`    any valid JuMP `@variable` syntax to define the value of the state variable at the beginning of the stage

# Examples
 * `@state(sp, 0 <= x[i=1:3] <= 1, x0==rand(3)[i] )`
 * `@state(sp,      y        <= 1, y0==0.5        )`
 * `@state(sp,      z            , z0==0.5        )`

"""
macro state(sp, x, x0)
    sp = esc(sp)                        # escape the model
    @assert x0.head == :call && x0.args[1] == :(==) # must be ==
    compsym, symin, rhs = x0.args                # name of the statein variable
    if is_comparison(x)
        if length(x.args) == 5          # double sided
            xin = _copy(x.args[3])       # variable is in middle
        elseif length(x.args) == 3      # single comparison
            xin = _copy(x.args[2])       # variable is second entry
        else
            error("Unknown format for $(x)")
        end
    else
        xin = _copy(x)  # no bounds
    end
    if isa(xin, Expr)                   # x has indices
        xin.args[1] = symin             # so just change the name
    else                                # its just a Symbol
        xin = symin                     # so change the Symbol
    end

    quote
        stateout = $(Expr(:macrocall, Symbol("@variable"), sp, esc(x)))
        statein  = $(Expr(:macrocall, Symbol("@variable"), sp, esc(xin), Expr(KW_SYM, START, esc(rhs))))
        statevariable!($sp, statein, stateout)
        stateout, statein
    end
end

"""
    @states(sp, begin
        stateleaving1, stateentering1
        stateleaving2, stateentering2
    end)

# Description

Define a new state variables in the subproblem `sp`.

# Arguments
* `sp`               the subproblem
* `stateleaving`     any valid JuMP `@variable` syntax to define the value of the state variable at the end of the stage
* `stateentering`    any valid JuMP `@variable` syntax to define the value of the state variable at the beginning of the stage

# Usage

    @states(sp, begin
        0 <= x[i=1:3] <= 1, x0==rand(3)[i]
             y        <= 1, y0==0.5
             z            , z0==0.5
     end)

"""
macro states(m, b)
    @assert b.head == :block || error("Invalid syntax for @states")
    code = quote end
    for line in b.args
        if !Base.Meta.isexpr(line, :line)
            if line.head == :tuple && length(line.args) == 2
                push!(code.args,
                    Expr(:macrocall, Symbol("@state"), esc(m), esc(line.args[1]), esc(line.args[2]))
                )
            else
                error("Unknown arguments in @states")
            end
        end
    end
    push!(code.args, :(nothing))
    return code
end
