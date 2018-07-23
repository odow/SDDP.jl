#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export Wasserstein

"""
    Wasserstein(alpha::Float64, solver::MathProgBase.AbstractMathProgSolver)

A distributionally-robust risk measure based on the Wasserstein distance.

As `alpha` increases, the measure becomes more risk-averse. When `alpha=0`, the
measure is equivalent to the expectation operator. As `alpha` increases, the
measure approaches the Worst-case risk measure.

This requires a `solver` that supports quadratic constraints.
"""
struct Wasserstein <: AbstractRiskMeasure
    alpha::Float64
    solver
    function Wasserstein(alpha::Float64, solver)
        if alpha < 0.0
            error("alpha must be in the range [0, ∞).")
        end
        return new(alpha, solver)
    end
end

function modifyprobability!(measure::Wasserstein,
    riskadjusted_distribution,
    original_distribution::Vector{Float64},
    observations::Vector{Float64},
    m::SDDPModel,
    sp::JuMP.Model
    )
    N = length(observations)
    wasserstein = JuMP.Model(solver=measure.solver)
    @variable(wasserstein, x[1:N] >= 0)
    @constraint(wasserstein, sum(x) == 1)
    @constraint(wasserstein, 
        # TODO: this should be ||x - original_distribution||₂ <= measure.alpha
        # but Clp and Ipopt (our solvers for the tests) don't support it. Is it
        # really a big deal?
        sum((x - original_distribution).^2) <= measure.alpha^2
    )
    if getsense(sp) == :Min
        @objective(wasserstein, Max, dot(observations, x))
    else
        @objective(wasserstein, Min, dot(observations, x))
    end
    solve(wasserstein)
    copy!(riskadjusted_distribution, getvalue(x))
end
