# Copyright 2017, Oscar Dowson

struct DematosCutOracle <: AbstractCutOracle
    sense::Symbol
    cuts::Vector{Cut}
    non_dominted_count::Vector{Int}
    statesvisited::Vector{Vector{Float64}}
    best_objectives::Vector{Float64}
    best_cut_idx::Vector{Float64}
end

DematosCutOracle(sense=:Min) = DematosCutOracle(sense, Cut[], Int[], Vector{Float64}[], Float64[], Float64[])

function storecut!(oracle::DematosCutOracle, m::SDDPModel, sp::JuMP.Model, cut::Cut)
    push!(oracle.cuts, cut)
    push!(oracle.non_dominated_count, 0)
    for (i, state) in enumerate(oracle.statesvisited)
        y = cut.intercept + dot(cut.coefficients, state)
        if (sense == :Min && y > best_objectives[i]) || (sense == :Max && y < best_objectives[i])
            best_objectives[i] = y
            non_dominated_count[best_cut_idx[i]] -= 1
            best_cut_idx[i] = length(oracle.cuts)
            non_dominated_count[best_cut_idx[i]] += 1
        end
    end

    current_state = stage(m, ext(sp).stage).state
    push!(oracle.statesvisited, copy(current_state))
    push!(oracle.best_objectives, cut.intercept + dot(cut.coefficients, current_state))
    push!(oracle.best_cut_idx, length(oracle.cuts))
    oracle.non_dominated_count[end] += 1
    for (i, cut) in enumerate(oracle.cuts)
        y = cut.intercept + dot(cut.coefficients, current_state)
        if (sense == :Min && y > best_objectives[end]) || (sense == :Max && y < best_objectives[end])
            best_objectives[end] = y
            non_dominated_count[best_cut_idx[end]] -= 1
            best_cut_idx[end] = i
            non_dominated_count[i] += 1
        end
    end
end

function validcuts(oracle::DematosCutOracle)
    active_cuts = Cut[]
    for cut in oracle.cuts
        if cut.non_dominated_count > 0
            push!(active_cuts, cut.cut)
        end
    end
    active_cuts
end
