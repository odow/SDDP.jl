struct Cut
    intercept::Float64
    duals::Vector{Float64}
end

# function constructcut!(
#     riskmeasure::AbstractRiskMeasure,
#     newprobability::Vector{Float64}     # probability of outcomes
#     oldprobability::Vector{Float64}     # probability of outcomes
#     m::SDDPModel,
#     stage::Int,
#     markovstate::Int,
#     x::Vector{Float64},             # value of state variables
#     pi::Vector{Vector{Float64}},    # vector of dual values
#     theta::Vector{Float64}          # vector of future value/cost
#     )
#
#
#     newprobabilities =
# end
