#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

function setnoise!(sp::JuMP.Model, noise::Noise)
    for (c, v) in zip(noise.constraints, noise.values)
        JuMP.setRHS(c, v)
    end
    if noise.has_objective_noise
        stageobjective!(valueoracle(sp), sp, noise.obj)
    end
end

function samplenoise(sp::JuMP.Model)
    noiseidx = sample(ext(sp).noiseprobability)
    return noiseidx, ext(sp).noises[noiseidx]
end

"""
    setnoiseprobability!(sp::JuMP.Model, distribution::Vector{Float64})

# Description

Set the probability distribution of the stagewise independent noise in the
`sp` subproblem.

# Arguments
* `sp`            the subproblem
* `distribution` vector containing the probability of each outcome occuring.
    Should sum to `1`. Defaults to the uniform distribution.

# Examples
If there are two realizations:
* `setnoiseprobability!(sp, [0.3, 0.7])`
* `setnoiseprobability!(sp, [0.5, 0.6])` will error!
"""
function setnoiseprobability!(sp::JuMP.Model, p::Vector{Float64})
    if abs(sum(p) - 1.0) > 1e-4 # check sum to one
        error("Discrete probability distribution of the noise must sum to one.")
    end
    resize!(ext(sp).noiseprobability, length(p))
    ext(sp).noiseprobability .= p
end

"""
    @rhsnoise(sp, rhs, constraint)

# Description

Add a constraint with a noise in the RHS vector to the subproblem `sp`.

# Arguments
* `sp`         the subproblem
* `rhs`        keyword argument `key=value` where `value` is a one-dimensional array containing the noise realisations
* `constraint` any valid JuMP `@constraint` syntax that includes the keyword defined by `rhs`

# Examples
* `@rhsnoise(sp, i=1:2, x + y <= i )`
* `@rhsnoise(sp, i=1:2, x + y <= 3 * rand(2)[i] )`
"""
macro rhsnoise(sp, kw, c)
    sp = esc(sp)                                # escape the model
    @assert kw.head == KW_SYM                   # check its a keyword
    noisevalues = esc(kw.args[2])            # get the vector of values
    @assert c.head == :call               # check c is a comparison constraint
    @assert length(c.args) == 3                 # check that it has (LHS, (comparison), RHS)
    @assert c.args[1]  in comparison_symbols # check valid constraint type
    constrexpr = :($(c.args[2]) - $(c.args[3])) # LHS - RHS
    quote
        rhs = Float64[]                         # intialise RHS vector
        $(esc(kw.args[1])) = $noisevalues[1] # initialise with first noise
        aff = copy($(esc(constrexpr)).coeffs)
        for val in $noisevalues    # for each noise
            $(esc(kw.args[1])) = val  # set the noisevalue
            if $(esc(constrexpr)).coeffs != aff
                error("""Uh oh! It looks like you've included the noise term in
                    constraint matrix. Unfortunately, the noise in @rhsnoise
                    can only appear in the RHS of the constraint. Upcoming
                    versions of JuMP may let us revisit this restriction in the
                    future.

                    As a work around, have a look at the asset_management.jl
                    example to see how constraint coefficients can be randomised
                    using a markov process.""")
            end
            push!(rhs, -$(esc(constrexpr)).constant)
         end
        $(esc(kw.args[1])) = $noisevalues[1] # initialise with first noise
        con = $(Expr(                           # add the constraint
                :macrocall, Symbol("@constraint"),
                sp,                             # the subproblem
                esc(c)                          # the constraint expression
                ))
        registernoiseconstraint!($sp, con, rhs)
        con
    end
end

function registernoiseconstraint!(sp::JuMP.Model, con::LinearConstraint, rhs::Vector{Float64})
    if length(ext(sp).noises) == 0
        for r in rhs
            push!(ext(sp).noises, Noise(false, AffExpr(), [con], [r]))
        end
    else
        @assert length(ext(sp).noises) == length(rhs)
        for (i, r) in enumerate(rhs)
            push!(ext(sp).noises[i].constraints, con)
            push!(ext(sp).noises[i].values, r)
        end
    end
end

"""
    @rhsnoises(sp, rhs, begin
        constraint
    end)

# Description

The plural form of `@rhsnoise` similar to the JuMP macro `@constraints`.

# Arguments

See `@rhsnoise`.

# Examples
    @rhsnoises(sp, i=1:2, begin
        x + y <= i
        x + y <= 3 * rand(2)[i]
    end)
"""
macro rhsnoises(m, kw, blk)
    @assert blk.head == :block || error("Invalid syntax for @rhsnoises")
    code = quote end
    for line in blk.args
        if !Base.Meta.isexpr(line, :line)
            if line.head == :call && line.args[1] in comparison_symbols
                push!(code.args,
                    Expr(:macrocall, Symbol("@rhsnoise"), esc(m), esc(kw), esc(line))
                )
            else
                error("Unknown arguments in @rhsnoises")
            end
        end
    end
    push!(code.args, :(nothing))
    return code
end

"""
    @stageobjective!(sp, kw=noises, objective)

# Description

Define an objective that depends on the realization of the stagewise noise.
`objective` can be any valid third argument to the JuMP `@objective` macro (i.e.
`@objective(sp, Min, objective)`) that utilises the variable `kw` that takes the
realizations defined in `noises`.

# Examples

    @stageobjective(sp, w=1:2, w * x)
    @stageobjective(sp, i=1:2, w[i]^2 * x)
    @stageobjective(sp, i=1:2, x[i])

"""
macro stageobjective(sp, kw, obj)
    sp = esc(sp)                                # escape the model
    @assert kw.head == KW_SYM                   # check its a keyword
    noisevalues = esc(kw.args[2])            # get the vector of values
    quote
        for (i, val) in enumerate($noisevalues)    # for each noise
            $(esc(kw.args[1])) = val  # set the noisevalue
            registernoiseobjective!($sp, $(esc(obj)), i)
         end
    end
end

macro stageobjective(sp, obj)
    quote
        stageobjective!($(esc(sp)), $(esc(obj)))
    end
end

function registernoiseobjective!(sp::JuMP.Model, objective, idx::Int)
    # cases. noises exist. in which case go through
    if length(ext(sp).noises) >= idx
        ext(sp).noises[idx].has_objective_noise = true
        ext(sp).noises[idx].obj = objective
    elseif idx == length(ext(sp).noises) + 1
        push!(ext(sp).noises, Noise(true, objective, [], []))
        # check we haven't added any constraint noises yet
        if length(ext(sp).noises[1].constraints) > 0
            error("You must have the same number of noises in the objective function as the constraint RHS's.")
        end
    end
end

include("deprecate.jl")
@deprecate_macro noise rhsnoise
@deprecate_macro noises rhsnoises
