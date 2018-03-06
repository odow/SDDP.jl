#  This function is copied from
# https://github.com/JuliaOpt/JuMP.jl/blob/963886fa8a9b668630f0b3704ed22c58a5947ee8/src/deprecated.jl

#############################################################################
#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

macro deprecate_macro(oldmacro,newmacro)
    oldmac = Symbol(string("@",oldmacro))
    newmac = Symbol(string("@",newmacro))
    s = string(oldmac," is deprecated, use ", newmac, " instead.")
    depwarn = :(Base.depwarn($s,$(Base.Meta.quot(oldmac))))
    @eval macro $oldmacro(args...)
        return Expr(:block, $depwarn, Expr(:macrocall, $(Base.Meta.quot(newmac)), [esc(x) for x in args]...))
    end
    eval(Expr(:export,oldmac))
    return
end

@deprecate_macro noise rhsnoise
@deprecate_macro noises rhsnoises
@deprecate stageobjective! setstageobjective!
@deprecate init! initializevaluefunction
@deprecate NestedAVaR(;lambda=1.0,beta=1.0) EAVaR(;lambda=lambda,beta=beta)
