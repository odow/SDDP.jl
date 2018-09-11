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
    model::SDDPModel,
    subproblem::JuMP.Model
    )
    N = length(observations)
    wasserstein = JuMP.Model(solver=measure.solver)
    @variables(wasserstein, begin
        transport_matrix[1:N, 1:N] >= 0
        adjusted_distribution[1:N] >= 0
        absolute_value[1:N, 1:N]   >= 0
    end)
    @constraints(wasserstein, begin
        [j=1:N], sum(transport_matrix[:, j]) == original_distribution[j]
        [i=1:N], sum(transport_matrix[i, :]) == adjusted_distribution[i]
        sum(transport_matrix[i, j] * abs(observations[i] - observations[j])
            for i in 1:N, j in 1:N) <= measure.alpha
    end)
    if getsense(subproblem) == :Min
        @objective(wasserstein, Max, dot(observations, adjusted_distribution))
    else
        @objective(wasserstein, Min, dot(observations, adjusted_distribution))
    end
    solve(wasserstein)
    copy!(riskadjusted_distribution, getvalue(adjusted_distribution))
end
