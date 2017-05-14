# Copyright 2017, Oscar Dowson
export DematosCutOracle

immutable DematosCutOracle <: AbstractCutOracle
    cuts::Vector{Cut}
    non_dominated_count::Vector{Int}
    statesvisited::Vector{Vector{Float64}}
    best_objectives::Vector{Float64}
    best_cut_idx::Vector{Int}
end

DematosCutOracle() = DematosCutOracle(Cut[], Int[], Vector{Float64}[], Float64[], Float64[])

function storecut!(o::DematosCutOracle, m::SDDPModel, sp::JuMP.Model, cut::Cut)
    sense = getsense(sp)
    push!(o.cuts, cut)
    push!(o.non_dominated_count, 0)
    for (i, state) in enumerate(o.statesvisited)
        y = cut.intercept + dot(cut.coefficients, state)
        if dominates(sense, y, o.best_objectives[i])
            o.best_objectives[i] = y
            o.non_dominated_count[o.best_cut_idx[i]] -= 1
            o.best_cut_idx[i] = length(o.cuts)
            o.non_dominated_count[o.best_cut_idx[i]] += 1
        end
    end

    # get the last state
    current_state = getstage(m, ext(sp).stage).state
    # add to oracle, and assume last cut is the best
    push!(o.statesvisited, copy(current_state))
    push!(o.best_objectives, cut.intercept + dot(cut.coefficients, current_state))
    push!(o.best_cut_idx, length(o.cuts))
    o.non_dominated_count[end] += 1
    # for all the cuts
    for (i, cut) in enumerate(o.cuts)
        # evaluate
        y = cut.intercept + dot(cut.coefficients, current_state)
        if dominates(sense, y, o.best_objectives[end])
            o.best_objectives[end] = y
            o.non_dominated_count[o.best_cut_idx[end]] -= 1
            o.best_cut_idx[end] = i
            o.non_dominated_count[i] += 1
        end
    end
end



function validcuts(o::DematosCutOracle)
    active_cuts = Cut[]
    for (cut, count) in zip(o.cuts, o.non_dominated_count)
        if count > 0
            push!(active_cuts, cut)
        end
    end
    active_cuts
end
