#  Copyright 2017-21, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This file implements two solver-specific functions need by BiObjectiveSDDP
# to compute the reduced costs of arbitrary objective vectors.
#
# See `BiObjectiveSDDP.get_BinvA` and `BiObjectiveSDDP.get_basis`.

import Gurobi
import SparseArrays

####
####    New functions for the Gurobi C API.
####

function get_basis(model::Gurobi.Optimizer)
    p = Ref{Cint}()
    @assert Gurobi.GRBgetintattr(model, "NumConstrs", p) == 0
    bhead = zeros(Cint, p[])
    ret = Gurobi.GRBgetBasisHead(model, bhead)
    @assert ret == 0
    bhead .+= 1
    return bhead
end

mutable struct GRBsvec
    len::Cint
    ind::Ptr{Cint}
    val::Ptr{Cdouble}
end

function get_BinvA(model::Gurobi.Optimizer)
    p = Ref{Cint}()
    @assert Gurobi.GRBgetintattr(model, "NumConstrs", p) == 0
    ncon = p[]
    @assert Gurobi.GRBgetintattr(model, "NumVars", p) == 0
    nvar = p[]
    function _GRBBinvRowi(model::Gurobi.Optimizer, i::Int)
        ind = zeros(Cint, ncon + nvar)
        val = zeros(Cdouble, ncon + nvar)
        x = GRBsvec(0, pointer(ind), pointer(val))
        GC.@preserve ind val x begin
            @assert Gurobi.GRBBinvRowi(model, i, pointer_from_objref(x)) == 0
        end
        return ind[1:x.len], val[1:x.len]
    end
    rows, cols, coefs = Cint[], Cint[], Cdouble[]
    for i in 1:ncon
        ind, val = _GRBBinvRowi(model, i - 1)
        append!(rows, fill(Cint(i), length(ind)))
        append!(cols, ind .+ 1)
        append!(coefs, val)
    end
    return SparseArrays.sparse(rows, cols, coefs)
end
