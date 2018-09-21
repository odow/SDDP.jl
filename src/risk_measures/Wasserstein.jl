#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export Wasserstein

"""
    Wasserstein(alpha::Float64, solver::MathProgBase.AbstractMathProgSolver)

An (incorrect) distributionally-robust risk measure based on the Wasserstein
distance.

As `alpha` increases, the measure becomes more risk-averse. When `alpha=0`, the
measure is equivalent to the expectation operator. As `alpha` increases, the
measure approaches the Worst-case risk measure.
"""
struct Wasserstein <: AbstractRiskMeasure
    alpha::Float64
    solver
    function Wasserstein(alpha::Float64, solver)
        if alpha < 0.0
            error("alpha cannot be $(alpha) as it must be in the range [0, ∞).")
        end
        return new(alpha, solver)
    end
end

function modify_probability(measure::Wasserstein,
                            riskadjusted_distribution,
                            original_distribution::Vector{Float64},
                            observations::Vector{Float64},
                            model::SDDPModel,
                            subproblem::JuMP.Model)
    # TODO(odow): this formulation is incorrect. Instead of using the
    # observations in the transport matrix, we  need to use the actual noise
    # values :'(.
    N = length(observations)
    wasserstein = JuMP.Model(solver=measure.solver)
    @variable(wasserstein, transport[1:N, 1:N] >= 0)
    @variable(wasserstein, adjusted_distribution[1:N] >= 0)
    for i in 1:N
        @constraint(wasserstein, sum(transport[:, i]) == original_distribution[i])
        @constraint(wasserstein, sum(transport[i, :]) == adjusted_distribution[i])
    end
    @constraint(wasserstein, sum(transport[i, j] * abs(observations[i] -
        observations[j]) for i in 1:N, j in 1:N) <= measure.alpha)
    if getsense(subproblem) == :Min
        @objective(wasserstein, Max, dot(observations, adjusted_distribution))
    else
        @objective(wasserstein, Min, dot(observations, adjusted_distribution))
    end
    solve(wasserstein)
    copy!(riskadjusted_distribution, getvalue(adjusted_distribution))
    return
end
