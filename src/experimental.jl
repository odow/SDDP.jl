#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
An example. I can't test on Travis as we can't install Gurobi.

    using SDDP, JuMP, Gurobi, Base.Test
    m = SDDPModel(
            sense           = :Min,
            stages          = 2,
            objective_bound = 0.0,
            solver          = GurobiSolver(OutputFlag=0) ) do sp, t
        @state(sp, x>=0, x0==1.0)
        r_x = addconstraintnoise!(sp, x0, [1,2])
        @constraint(sp, x == r_x)
        @stageobjective(sp, x)
    end
    solve(m, iteration_limit=5)
    @test getbound(m) == 3.75

=#

"""
    addconstraintnoise!(sp::JuMP.Model, x::JuMP.Variable, realizations::Vector{T}) where T

Add a stagewise-independent coefficient term into the model.

Important:
  - Requires a solver that has implemented the `MathProgBase.changecoeffs!`
    method.
  - Cannot be used if `JuMP.solvehook` is set (e.g., by SDDiP.jl).

If multiple terms have random coefficents, then the length of realizations must
be the same in every case. This can also be used with stagewise-independent
right-hand side terms (via `@rhsnoise`), again with equal numbers of
realizations.

### Example
If `w ∈ [1,2,3]` with equal probability, and we want to add the constraint:
    @constraint(sp, 2x + w*x <= 1)
Then
    wx = addconstraintnoise!(sp, x, [1,2,3])
    @constraint(m, 2x + wx <= 1)

Note: This is a bit of a hack. Given a term `w*x` in the model, for example, in
a constraint such as `2x + w*x <= 1`, we introduce a dummy variable and a dummy
constraint so that: `dummyconstraint: dummy_variable == w * x`. This allows us
to deterministically work around the fact that JuMP cannot modify the
coefficient matrix.
"""
function addconstraintnoise!(sp, x, realizations)
    if !haskey(sp.ext, :changes)
        # initialize storage
        sp.ext[:changes] = (Int[], Int[], Int[])
    end
    # create an anonymous variable that is `noise * x`
    noise_state = @variable(sp)
    # create a nicely canonicalized constraint `noise_state == noise * x` with
    # the dummy value of `noise=1.0` to start
    c = @constraint(sp, 1.0 * x - noise_state == 0.0)
    # create a dummy variable for the noise
    w = @variable(sp)
    # leverage SDDP.jl machinery to perform the random variable stuff,
    # substituting in the dummy variable for a fake constraint.
    c2 = @rhsnoise(sp, i=realizations, w==i)
    # store information for solvehook
    push!(sp.ext[:changes][1],  c.idx)
    push!(sp.ext[:changes][2],  x.col)
    push!(sp.ext[:changes][3], c2.idx)
    # set solvehook if not already
    if sp.solvehook == nothing
        function our_solve_hook(sp; kwargs...)
            rows = sp.ext[:changes][1]::Vector{Int}
            cols = sp.ext[:changes][2]::Vector{Int}
            noise_terms = fill(0.0, length(rows))
            for (i, w) in enumerate(sp.ext[:changes][3])
                noise_terms[i] = sp.linconstr[w].ub
            end
            JuMP.MathProgBase.changecoeffs!(JuMP.internalmodel(sp),
                rows, cols, noise_terms)
            solve(sp, ignore_solve_hook=true)
        end
        JuMP.setsolvehook(sp, our_solve_hook)
        JuMP.build(sp)
    end
    # return the anonymous variable that is `noise * x`
    return noise_state
end
