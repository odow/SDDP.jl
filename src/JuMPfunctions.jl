#=
    This function is a modified version of that contained in the file
https://github.com/JuliaOpt/JuMP.jl/blob/223103c36e976ef670cce15591cb337d5d475a68/src/solvers.jl

    It was simplified for SDDP.jl use

=#
#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/solvers.jl
# Handles conversion of the JuMP Model into a format that can be passed
# through the MathProgBase interface to solvers, and ongoing updating of
# that representation if supported by the solver.
#############################################################################

function jumpsolve(m::JuMP.Model; suppress_warnings=false,
                ignore_solve_hook=(m.solvehook===nothing),
                relaxation=false,
                kwargs...)
    # If the user or an extension has provided a solve hook, call
    # that instead of solving the model ourselves
    if !ignore_solve_hook
        return m.solvehook(m; suppress_warnings=suppress_warnings, kwargs...)::Symbol
    end

    isempty(kwargs) || error("Unrecognized keyword arguments: $(join([k[1] for k in kwargs], ", "))")

    # Clear warning counters
    m.map_counter = 0
    m.operator_counter = 0

    # Analyze the problems traits to determine what solvers we can use
    traits = JuMP.ProblemTraits(m, relaxation=relaxation)

    # Build the MathProgBase model from the JuMP model
    JuMP.build(m, traits=traits, suppress_warnings=suppress_warnings, relaxation=relaxation)

    # Solve the problem
    JuMP.MathProgBase.optimize!(m.internalModel)
    stat::Symbol = JuMP.MathProgBase.status(m.internalModel)

    # Extract solution from the solver
    numRows, numCols = length(m.linconstr), m.numCols
    m.objBound = NaN
    m.objVal = NaN
    if length(m.colVal) != numCols
        m.colVal = fill(NaN, numCols)
    end
    if length(m.linconstrDuals) != numCols
        m.linconstrDuals = fill(NaN, numCols)
    end

    discrete = !relaxation && (traits.int || traits.sos)
    if stat == :Optimal && !discrete
        m.linconstrDuals = JuMP.MathProgBase.getconstrduals(m.internalModel)[1:numRows]
    else
        # Problem was not solved to optimality, attempt to extract useful
        # information anyway
        suppress_warnings || warn("Not solved to optimality, status: $stat")
    end

    # If the problem was solved, or if it terminated prematurely, try
    # to extract a solution anyway. This commonly occurs when a time
    # limit or tolerance is set (:UserLimit)
    if !(stat == :Infeasible || stat == :Unbounded)
        m.objVal = JuMP.MathProgBase.getobjval(m.internalModel) + m.obj.aff.constant
        m.colVal .= JuMP.MathProgBase.getsolution(m.internalModel)[1:numCols]
    end
    # don't keep relaxed model in memory
    relaxation && (m.internalModelLoaded = false)
    # Return the solve status
    stat
end
