#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

function setscenario!(sp::JuMP.Model, scenario::Scenario)
    for (c, v) in zip(scenario.constraints, scenario.values)
        JuMP.setRHS!(c, v)
    end
end

"""
    setscenarioprobability!(sp::JuMP.Model, p::Vector{Float64})
"""
function setscenarioprobability!(sp::JuMP.Model, p::Vector{Float64})
    @assert abs(sum(p) - 1.0) < 1e-4 # check sum to one
    resize!(ext(sp).scenarioprobability, length(p))
    ext(sp).scenarioprobability .= p
end

"""
    @scenario(sp, rhs, constraint)
Add a scenario constraint (changes in RHS vector) to the subproblem `sp`.
Arguments:
    sp             the subproblem
    rhs            keyword argument `key=value` where `value` is a one-dimensional array containing the scenario realisations
    constraint     any valid JuMP `@constraint` syntax that includes the keyword defined by `rhs`
Usage:
    @scenario(sp, i=1:2, x + y <= i )
    @scenario(sp, i=1:2, x + y <= 3 * rand(2)[i] )
"""
macro scenario(sp, kw, c)
    sp = esc(sp)                                # escape the model
    @assert kw.head == :(=)                     # check its a keyword
    scenariovalues = esc(kw.args[2])            # get the vector of values
    @assert c.head == :call               # check c is a comparison constraint
    @assert length(c.args) == 3                 # check that it has (LHS, (comparison), RHS)
    @assert c.args[1]  in comparison_symbols # check valid constraint type
    constrexpr = :($(c.args[2]) - $(c.args[3])) # LHS - RHS
    quote
        rhs = Float64[]                         # intialise RHS vector
        for scenariovalue in $scenariovalues    # for each scenario
            $(esc(kw.args[1])) = scenariovalue  # set the scenariovalue
            push!(rhs, -$(esc(constrexpr)).constant)
         end
        $(esc(kw.args[1])) = $scenariovalues[1] # initialise with first scenario
        con = $(Expr(                           # add the constraint
                :macrocall, Symbol("@constraint"),
                sp,                             # the subproblem
                esc(c)                          # the constraint expression
                ))
        registerscenarioconstraint!($sp, con, rhs)
        con
    end
end

function registerscenarioconstraint!(sp::JuMP.Model, con::LinearConstraint, rhs::Vector{Float64})
    if length(ext(sp).scenarios) == 0
        for r in rhs
            push!(ext(sp).scenarios, Scenario([con], [r]))
        end
    else
        @assert length(ext(sp).scenarios) == length(rhs)
        for (i, r) in enumerate(rhs)
            push!(ext(sp).scenarios[i].constraints, con)
            push!(ext(sp).scenarios[i].values, r)
        end
    end
end

"""
    @scenarios(sp, rhs, begin
        constraint
    end)
The plural form of `@scenario` similar to the JuMP macro `@constraints`.
Usage:
    @scenario(sp, i=1:2, begin
               x + y <= i
               x + y <= 3 * rand(2)[i]
    end)
"""
macro scenarios(m, kw, blk)
    @assert blk.head == :block || error("Invalid syntax for @scenario")
    code = quote end
    for line in blk.args
        if !Base.Meta.isexpr(line, :line)
            if line.head == :call && line.args[1] in comparison_symbols
                push!(code.args,
                    Expr(:macrocall, Symbol("@scenario"), esc(m), esc(kw), esc(line))
                )
            else
                error("Unknown arguments in @scenario")
            end
        end
    end
    push!(code.args, :(nothing))
    return code
end
