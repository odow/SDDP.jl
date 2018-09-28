abstract type AbstractRiskMeasure end

"""
    adjust_probability(measure::Expectation
                       risk_adjusted_probability::Vector{Float64},
                       original_probability::Vector{Float64},
                       noise_support::Vector{T},
                       objective_realizations::Vector{Float64},
                       is_minimization::Bool)
"""
function adjust_probability end

# ========================== The Expectation Operator ==========================

struct Expectation <: AbstractRiskMeasure end

function adjust_probability(measure::Expectation,
                            risk_adjusted_probability::Vector{Float64},
                            original_probability::Vector{Float64},
                            noise_support::Vector{T},
                            objective_realizations::Vector{Float64},
                            is_minimization::Bool) where T
    risk_adjusted_probability .= original_probability
    return
end

# ========================== The Worst Case Operator ===========================

struct WorstCase <: AbstractRiskMeasure end

function adjust_probability(measure::WorstCase,
                            risk_adjusted_probability::Vector{Float64},
                            original_probability::Vector{Float64},
                            noise_support::Vector{T},
                            objective_realizations::Vector{Float64},
                            is_minimization::Bool) where T
    risk_adjusted_probability .= 0.0
    worst_index = 1
    worst_observation = is_minimization ? -Inf : Inf
    for (index, (probability, observation)) in enumerate(
            zip(original_probability, objective_realizations))
        if probability > 0.0
            if (is_minimization && observation > worst_observation) ||
                    (!is_minimization && observation < worst_observation)
                worst_index = index
                worst_observation = observation
            end
        end
    end
    risk_adjusted_probability[worst_index] = 1.0
    return
end
