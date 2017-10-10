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
    @timeit TIMER "prep JuMP model" begin
        # Clear warning counters
        m.map_counter = 0
        m.operator_counter = 0

        # Remember if the solver was initially unset so we can restore
        # it to be unset later
        unset = m.solver == JuMP.UnsetSolver()

        # Analyze the problems traits to determine what solvers we can use
        traits = JuMP.ProblemTraits(m, relaxation=relaxation)

        # Build the MathProgBase model from the JuMP model
        JuMP.build(m, traits=traits, suppress_warnings=suppress_warnings, relaxation=relaxation)

        # If the model is a general nonlinear, use different logic in
        # nlp.jl to solve the problem
        traits.nlp && return JuMP.solvenlp(m, traits, suppress_warnings=suppress_warnings)
    end
    # Solve the problem
    @timeit TIMER "optimize!" begin
        JuMP.MathProgBase.optimize!(m.internalModel)
    end

    @timeit TIMER "getsolution" begin
        stat::Symbol = JuMP.MathProgBase.status(m.internalModel)

        # Extract solution from the solver
        numRows, numCols = length(m.linconstr), m.numCols
        m.objBound = NaN
        m.objVal = NaN
        m.colVal = fill(NaN, numCols)
        m.linconstrDuals = Array{Float64}(0)

        discrete = !relaxation && (traits.int || traits.sos)
        if stat == :Optimal
            # If we think dual information might be available, try to get it
            # If not, return an array of the correct length
            if discrete
                m.redCosts = fill(NaN, numCols)
                m.linconstrDuals = fill(NaN, numRows)
            else
                if !traits.conic
                    m.redCosts = try
                        JuMP.MathProgBase.getreducedcosts(m.internalModel)[1:numCols]
                    catch
                        fill(NaN, numCols)
                    end

                    m.linconstrDuals = try
                        JuMP.MathProgBase.getconstrduals(m.internalModel)[1:numRows]
                    catch
                        fill(NaN, numRows)
                    end
                elseif !traits.qp && !traits.qc
                    fillConicDuals(m)
                end
            end
        else
            # Problem was not solved to optimality, attempt to extract useful
            # information anyway
            suppress_warnings || warn("Not solved to optimality, status: $stat")
            # Some solvers provide infeasibility rays (dual) or unbounded
            # rays (primal) for linear problems. Store these as the solution
            # if the exist.
            if traits.lin
                if stat == :Infeasible
                    m.linconstrDuals = try
                        infray = JuMP.MathProgBase.getinfeasibilityray(m.internalModel)
                        @assert length(infray) == numRows
                        infray
                    catch
                        suppress_warnings || warn("Infeasibility ray (Farkas proof) not available")
                        fill(NaN, numRows)
                    end
                elseif stat == :Unbounded
                    m.colVal = try
                        unbdray = JuMP.MathProgBase.getunboundedray(m.internalModel)
                        @assert length(unbdray) == numCols
                        unbdray
                    catch
                        suppress_warnings || warn("Unbounded ray not available")
                        fill(NaN, numCols)
                    end
                end
            end
            # conic duals (currently, SOC and SDP only)
            if !discrete && traits.conic && !traits.qp && !traits.qc
                if stat == :Infeasible
                    fillConicDuals(m)
                end
            end
        end

        # If the problem was solved, or if it terminated prematurely, try
        # to extract a solution anyway. This commonly occurs when a time
        # limit or tolerance is set (:UserLimit)
        if !(stat == :Infeasible || stat == :Unbounded)
            if discrete
                try
                    # Do a separate try since getobjval could work while getobjbound does not and vice versa
                    objBound = JuMP.MathProgBase.getobjbound(m.internalModel) + m.obj.aff.constant
                    m.objBound = objBound
                end
            end
            try
                objVal = JuMP.MathProgBase.getobjval(m.internalModel) + m.obj.aff.constant
                colVal = JuMP.MathProgBase.getsolution(m.internalModel)[1:numCols]
                # Rescale off-diagonal terms of SDP variables
                if traits.sdp
                    offdiagvars = offdiagsdpvars(m)
                    colVal[offdiagvars] /= sqrt(2)
                end
                # Don't corrupt the answers if one of the above two calls fails
                m.objVal = objVal
                m.colVal = colVal
            end
        end

        # The MathProgBase interface defines a conic problem to always be
        # a minimization problem, so we need to flip the objective before
        # reporting it to the user
        if traits.conic && m.objSense == :Max
            m.objBound *= -1
            m.objVal *= -1
        end

        # If the solver was initially not set, we will restore this status
        # and drop the internal MPB model. This is important for the case
        # where the solver used changes between solves because the user
        # has changed the problem class (e.g. LP to MILP)
        if unset
            m.solver = JuMP.UnsetSolver()
            if traits.int
                m.internalModelLoaded = false
            end
        end

        # don't keep relaxed model in memory
        relaxation && (m.internalModelLoaded = false)
    end
    # Return the solve status
    stat
end
