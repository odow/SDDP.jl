#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

# ========================== The Expectation Operator ======================== #

"""
    Expectation()

The Expectation risk measure. Identical to taking the expectation with respect
to the nominal distribution.
"""
struct Expectation <: AbstractRiskMeasure end

function adjust_probability(
    measure::Expectation,
    risk_adjusted_probability::Vector{Float64},
    original_probability::Vector{Float64},
    noise_support::Vector,
    objective_realizations::Vector{Float64},
    is_minimization::Bool,
)
    risk_adjusted_probability .= original_probability
    return 0.0
end

# ========================== The Worst Case Operator ========================= #

"""
    WorstCase()

The worst-case risk measure. Places all of the probability weight on the worst
outcome.
"""
struct WorstCase <: AbstractRiskMeasure end

function adjust_probability(
    measure::WorstCase,
    risk_adjusted_probability::Vector{Float64},
    original_probability::Vector{Float64},
    noise_support::Vector,
    objective_realizations::Vector{Float64},
    is_minimization::Bool,
)
    risk_adjusted_probability .= 0.0
    worst_index = 1
    worst_observation = is_minimization ? -Inf : Inf
    for (index, (probability, observation)) in
        enumerate(zip(original_probability, objective_realizations))
        if probability > 0.0
            if (is_minimization && observation > worst_observation) ||
               (!is_minimization && observation < worst_observation)
                worst_index = index
                worst_observation = observation
            end
        end
    end
    risk_adjusted_probability[worst_index] = 1.0
    return 0.0
end

# =================================== AV@R =================================== #

"""
    AVaR(β)

The average value at risk (AV@R) risk measure.

Computes the expectation of the β fraction of worst outcomes. β must be in `[0,
1]`. When `β=1`, this is equivalent to the [`Expectation`](@ref) risk measure.
When `β=0`, this is equivalent  to the [`WorstCase`](@ref) risk measure.

AV@R is also known as the conditional value at risk (CV@R) or expected
shortfall.
"""
struct AVaR <: AbstractRiskMeasure
    β::Float64
    function AVaR(β::Float64)
        if !(0 <= β <= 1)
            throw(ArgumentError("Risk-quantile β must be in [0, 1]. Currently it is $(β)."))
        end
        return new(β)
    end
end

function adjust_probability(
    measure::AVaR,
    risk_adjusted_probability::Vector{Float64},
    original_probability::Vector{Float64},
    noise_support::Vector,
    objective_realizations::Vector{Float64},
    is_minimization::Bool,
)
    if measure.β ≈ 0.0
        return adjust_probability(
            WorstCase(),
            risk_adjusted_probability,
            original_probability,
            noise_support,
            objective_realizations,
            is_minimization,
        )
    elseif measure.β ≈ 1.0
        return adjust_probability(
            Expectation(),
            risk_adjusted_probability,
            original_probability,
            noise_support,
            objective_realizations,
            is_minimization,
        )
    end
    risk_adjusted_probability .= 0.0
    quantile_collected = 0.0
    for i in sortperm(objective_realizations, rev = is_minimization)
        quantile_collected >= measure.β && break
        avar_prob = min(original_probability[i], measure.β - quantile_collected) / measure.β
        risk_adjusted_probability[i] = avar_prob
        quantile_collected += avar_prob * measure.β
    end
    return 0.0
end

# ============================ ConvexCombination ============================= #

"""
    ConvexCombination((weight::Float64, measure::AbstractRiskMeasure)...)

Create a weighted combination of risk measures.

### Examples

    SDDP.ConvexCombination(
        (0.5, SDDP.Expectation()),
        (0.5, SDDP.AVaR(0.25))
    )

Convex combinations can also be constructed by adding weighted risk measures
together as follows:

    0.5 * SDDP.Expectation() + 0.5 * SDDP.AVaR(0.5)
"""
struct ConvexCombination{T<:Tuple} <: AbstractRiskMeasure
    measures::T
end
ConvexCombination(args::Tuple...) = ConvexCombination(args)
function Base.show(io::IO, measure::ConvexCombination)
    print(io, "A convex combination of ")
    is_first = true
    for m in measure.measures
        !is_first && print(io, " + ")
        print(io, m[1], " * ", m[2])
        is_first = false
    end
end
import Base: +, *

function Base.:+(a::ConvexCombination, b::ConvexCombination)
    return ConvexCombination(a.measures..., b.measures...)
end

function Base.:*(lhs::Float64, rhs::AbstractRiskMeasure)
    return ConvexCombination(((lhs, rhs),))
end

function adjust_probability(
    measure::ConvexCombination,
    risk_adjusted_probability::Vector{Float64},
    original_probability::Vector{Float64},
    noise_support::Vector,
    objective_realizations::Vector{Float64},
    is_minimization::Bool,
)
    risk_adjusted_probability .= 0.0
    partial_distribution = similar(risk_adjusted_probability)
    for (weight, measure) in measure.measures
        partial_distribution .= 0.0
        adjust_probability(
            measure,
            partial_distribution,
            original_probability,
            noise_support,
            objective_realizations,
            is_minimization,
        )
        risk_adjusted_probability .+= weight * partial_distribution
    end
    return 0.0
end

# =================================== EAV@R ================================== #

"""
    EAVaR(;lambda=1.0, beta=1.0)

A risk measure that is a convex combination of Expectation and Average Value @
Risk (also called Conditional Value @ Risk).

        λ * E[x] + (1 - λ) * AV@R(1-β)[x]

### Keyword Arguments

* `lambda`: Convex weight on the expectation (`(1-lambda)` weight is put on the
  AV@R component. Inreasing values of `lambda` are less risk averse (more
  weight on expectation).

* `beta`: The quantile at which to calculate the Average Value @ Risk.
  Increasing values of `beta` are less risk averse. If `beta=0`, then the AV@R
  component is the worst case risk measure.
"""
function EAVaR(; lambda::Float64 = 1.0, beta::Float64 = 1.0)
    if !(0.0 <= lambda <= 1.0)
        error(
            "Lambda must be in the range [0, 1]. Increasing values of " *
            "lambda are less risk averse. lambda=1 is identical to taking " *
            "the expectation.",
        )
    end
    if !(0.0 <= beta <= 1.0)
        error(
            "Beta must be in the range [0, 1]. Increasing values of beta " *
            "are less risk averse. lambda=1 is identical to taking the " *
            "expectation.",
        )
    end
    return lambda * Expectation() + (1 - lambda) * AVaR(beta)
end

# ================================= Modified-Χ² ============================== #

#=
This code was contributed by Lea Kapelevich.

In a Distributionally Robust Optimization (DRO) approach, we modify the
probabilities we associate with all future scenarios so that the resulting
probability distribution is the "worst case" probability distribution, in some
sense.

In each backward pass we will compute a worst case probability distribution
vector ̃p. We compute ̃p so that:

̄p ∈ argmax{̃pᵀ̃z}
    ||̃p - ̃a||₂ ≤ r
    ∑̃p = 1
    ̃p ≥ 0

where

 1. ̃z is a vector of future costs. We assume that our aim is to minimize
    future cost pᵀ̃z. If we maximize reward, we would have ̃p ∈ argmin{̃pᵀ̃z}.
2. ̄a is the uniform distribution
3. r is a user specified radius - the larger the radius, the more conservative
   the policy.

Note: the largest radius that will work with S scenarios is sqrt((S-1)/S).
=#

"""
    ModifiedChiSquared(radius::Float64; minimum_std=1e-5)

The distributionally robust SDDP risk measure of
Philpott, A., de Matos, V., Kapelevich, L. Distributionally robust SDDP.
Computational Management Science (2018) 165:431-454.

If the uncorrected standard deviation of the objecive realizations is less than
`minimum_std`, then the risk-measure will default to `Expectation()`.
"""
struct ModifiedChiSquared <: AbstractRiskMeasure
    radius::Float64
    minimum_std::Float64
    function ModifiedChiSquared(radius::Float64; minimum_std::Float64 = 1e-5)
        if abs(radius) < 1e-9
            @warn(
                "Radius is very small. You should probably use " *
                "`SDDP.Expectation()` instead."
            )
        end
        return new(radius, minimum_std)
    end
end

function Base.show(io::IO, measure::ModifiedChiSquared)
    print(io, "ModifiedChiSquared with radius=$(measure.radius)")
end

function adjust_probability(
    measure::ModifiedChiSquared,
    risk_adjusted_probability::Vector{Float64},
    original_probability::Vector{Float64},
    noise_support::Vector,
    objective_realizations::Vector{Float64},
    is_minimization::Bool,
)
    if Statistics.std(objective_realizations, corrected = false) < measure.minimum_std
        return adjust_probability(
            Expectation(),
            risk_adjusted_probability,
            original_probability,
            noise_support,
            objective_realizations,
            is_minimization,
        )
    end
    m = length(objective_realizations)
    if all(x -> x ≈ 1 / m, original_probability)
        uniform_dro(
            measure,
            risk_adjusted_probability,
            original_probability,
            objective_realizations,
            is_minimization,
        )
    else
        non_uniform_dro(
            measure,
            risk_adjusted_probability,
            original_probability,
            objective_realizations,
            is_minimization,
        )
    end
    return 0.0
end

"""
Algorithm (1) of Philpott et al. Assumes that the nominal distribution is _not_
uniform.
"""
function non_uniform_dro(
    measure::ModifiedChiSquared,
    p::Vector{Float64},
    q::Vector{Float64},
    z::Vector{Float64},
    is_minimization::Bool,
)
    m = length(z)
    if Statistics.std(z) < 1e-6
        p .= q
        return 0.0
    end
    if !is_minimization
        z = -z
    end
    # step 1
    K = collect(1:m)
    # check if nomial probability is 0
    for i in K
        if isapprox(q[i], 0.0, atol = 1e-10)
            p[i] = 0
            splice!(K, i)
        end
    end
    #update m
    m = length(K)
    # use this to store the index popped out of K
    not_in_K = Int[]
    # step 2
    while length(K) > 1
        # step 2(a)
        z_bar = sum(z[i] for i in K) / length(K)
        s = sqrt(sum(z[i]^2 - z_bar^2 for i in K) / length(K))
        if isapprox(s, 0.0, atol = 1e-10)
            error("s is too small")
	end
        # step 2(b)
        if length(K) == m
            for i in K
                p[i] = q[i] + (z[i] - z_bar) / (sqrt(m) * s) * measure.radius
            end
        else
            for i in not_in_K
                p[i] = 0
            end

            sum_qj = sum(q[i] for i in not_in_K)
            sum_qj_squared = sum(q[i]^2 for i in not_in_K)
            len_k = length(K)
            n = sqrt(len_k * (measure.radius^2 - sum_qj_squared)-sum_qj^2)
            for i in K
                p[i] = q[i] + 1 / len_k * (sum_qj + n * (z[i] - z_bar) / s)
            end
        end

        # step 2(c)
        if all(pi -> pi >= 0.0, p)
            return 0.0
        end

        # find i(K)
        # find the list of indexes for which p is less than 0
        negative_p = K[ p[K] .< 0]
        computed_r = zeros(0)
        sum_qj = 0
        sum_qj_squared = 0
        if length(not_in_K) > 0
            sum_qj = sum(q[i] for i in not_in_K)
            sum_qj_squared = sum(q[i]^2 for i in not_in_K)
        end
        len_k = length(K)
        computed_r = [
            (((-q[i] * len_k  - sum_qj)/((z[i] - z_bar))/s)^2 + sum_qj_squared^2)/len_k + sum_qj_squared
            for i in negative_p
        ]
        i_K = negative_p[argmin(computed_r)]
        append!(not_in_K, i_K)
        filter!(e -> e != i_K,K)
    end
    # step 3
    for i in not_in_K
        p[i] = 0
    end
    p[K[1]] = 1
    return 0.0
end


"""
Algorithm (2) of Philpott et al. Assumes that the nominal distribution is
uniform.
"""
function uniform_dro(
    measure::ModifiedChiSquared,
    risk_adjusted_probability::Vector{Float64},
    original_probability::Vector{Float64},
    objective_realizations::Vector{Float64},
    is_minimization::Bool,
)
    m = length(objective_realizations)
    # Take a permuted view of `risk_adjusted_probability` so we can refer to
    # `p[i]` instead of `risk_adjusted_probability[perm[i]]`.
    perm = sortperm(objective_realizations, rev = !is_minimization)
    p = view(risk_adjusted_probability, perm)
    z = view(objective_realizations, perm)
    # Compute the new probabilities according to Algorithm (2) of the Philpott
    # et al. paper.
    # Step (1):
    @inbounds for k = 0:m-2
        # Step (1a):
        z_bar = sum(z[i] for i = (k+1):m) / (m - k)
        s² = sum(z[i]^2 - z_bar^2 for i = (k+1):m) / (m - k)
        # Due to numerical error, s² may sometimes be a little bit negative.
        if s² < -1e-8
            error("Something unexpected happened with s² term: `$(s²) < 0.0`.")
        elseif s² <= 0.0
            error("`s²<0`: choose a larger threshold for `minimum_std`.")
        end
        # Step (1b): note that we cache a couple of terms that don't depend on i
        #            to speed things up.
        term_1 = 1 / (m - k)
        term_2 = sqrt((m - k) * measure.radius^2 - k / m) / ((m - k) * sqrt(s²))
        # We really should set p[i] = 0 for i = 1, ..., k. But since we don't
        # touch p[k-1] again, we can just set the k'th element to 0.
        if k > 0
            p[k] = 0.0
        end
        if is_minimization
            @inbounds for i = (k+1):m
                p[i] = term_1 + term_2 * (z[i] - z_bar)
            end
        else
            # Okay, here's the rub: we should have converted
            # objective_realizations (rewards) into costs by negating them. This
            # would have required a copy. This means that z_bar is in fact the
            # -ve of what it should be. `s` is fine since it is a difference of
            # squares. Thus, all we have to do is negate both z[i] and z_bar
            # here.
            @inbounds for i = (k+1):m
                p[i] = term_1 + term_2 * (z_bar - z[i])
            end
        end
        # Step (1c)
        if p[k+1] >= 0.0
            return 0.0
        end
    end
    # Step (2):
    p[end] = 1.0
    return 0.0
end

# ================================= Wasserstein ============================== #

"""
    Wasserstein(norm::Function, solver_factory; alpha::Float64)

A distributionally-robust risk measure based on the Wasserstein distance.

As `alpha` increases, the measure becomes more risk-averse. When `alpha=0`, the
measure is equivalent to the expectation operator. As `alpha` increases, the
measure approaches the Worst-case risk measure.
"""
struct Wasserstein{T,F} <: AbstractRiskMeasure
    alpha::Float64
    solver_factory::T
    norm::F
    function Wasserstein(norm::Function, solver_factory; alpha::Float64)
        if alpha < 0.0
            error("alpha cannot be $(alpha) as it must be in the range [0, ∞).")
        end
        return new{typeof(solver_factory),typeof(norm)}(alpha, solver_factory, norm)
    end
end
Base.show(io::IO, measure::Wasserstein) = print(io, "SDDP.Wasserstein")

function adjust_probability(
    measure::Wasserstein,
    risk_adjusted_probability::Vector{Float64},
    original_probability::Vector{Float64},
    noise_support::Vector,
    objective_realizations::Vector{Float64},
    is_minimization::Bool,
)
    N = length(objective_realizations)
    wasserstein = JuMP.Model(measure.solver_factory)
    @variable(wasserstein, z[1:N, 1:N] >= 0)
    @variable(wasserstein, p[1:N] >= 0)
    for i = 1:N
        @constraint(wasserstein, sum(z[:, i]) == original_probability[i])
        @constraint(wasserstein, sum(z[i, :]) == p[i])
    end
    @constraint(
        wasserstein,
        sum(
            measure.norm(noise_support[i], noise_support[j]) * z[i, j] for i = 1:N, j = 1:N
        ) <= measure.alpha
    )
    objective_sense = is_minimization ? MOI.MAX_SENSE : MOI.MIN_SENSE
    @objective(
        wasserstein,
        objective_sense,
        sum(objective_realizations[i] * p[i] for i = 1:N)
    )
    JuMP.optimize!(wasserstein)
    if JuMP.primal_status(wasserstein) != MOI.FEASIBLE_POINT
        error(
            "Unable to solver Wasserstein subproblem. Status: ",
            JuMP.termination_status(wassserstein),
        )
    end
    copyto!(risk_adjusted_probability, JuMP.value.(p))
    return 0.0
end

# ================================= Entropic ============================== #

"""
    Entropic(γ::Float64)

The entropic risk measure as described by Dowson, Morton, and Pagnoncelli
(2020). Multistage stochastic programs with the entropic risk measure.
http://www.optimization-online.org/DB_HTML/2020/08/7984.html
"""
mutable struct Entropic <: AbstractRiskMeasure
    γ::Float64
end

function Base.show(io::IO, ent::Entropic)
    return print(io, "Entropic risk measure with γ = $(ent.γ)")
end

function adjust_probability(
    measure::Entropic,
    Q::Vector{Float64},
    p::Vector{Float64},
    ::Vector,
    X::Vector{Float64},
    is_min::Bool,
)
    if measure.γ == 0.0  # Special case of entropic: if γ = 0, F[X] = E[X].
        Q .= p
        return 0.0
    end
    # Handle maximization problems by negating γ. Usually we would negate X, but
    # with convex RM's, we also have to negate the α(q) term, so it's easier to
    # just negate γ.
    γ = (is_min ? 1.0 : -1.0) * measure.γ
    # Use `BigFloat` to avoid overflow that occurs when calculating `exp(x)`.
    y = p .* exp.(big.(γ .* X))
    Q .= y / sum(y)
    α = sum(
        # To avoid numerical issues with the log, skip elements that are `≈ 0`.
        qi * log(qi / pi) for (pi, qi) in zip(p, Q) if pi > 1e-14  && qi > 1e-14
    )
    return -α / γ
end
