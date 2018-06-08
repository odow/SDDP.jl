#  Copyright 2017, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

mutable struct SampledState
    state::Vector{Float64}
    best_objective::Float64
    best_cut_index::Int
end

mutable struct StoredCut
    cut::Cut
    non_dominated_count::Int
end

"""
    LevelOneCutOracle()

# Description

Initialize the cut oracle for Level One cut selection. See:

V. de Matos, A. Philpott, E. Finardi, Improving the performance of Stochastic
Dual Dynamic Programming, Journal of Computational and Applied Mathematics
290 (2015) 196–208.
"""
mutable struct LevelOneCutOracle <: AbstractCutOracle
    cuts::Vector{StoredCut}
    states::Vector{SampledState}
    sampled_states::Set{Vector{Float64}}
    LevelOneCutOracle() = new(StoredCut[], SampledState[], Set{Vector{Float64}}())
end

DematosCutOracle() = LevelOneCutOracle()

function storecut!(o::LevelOneCutOracle, m::SDDPModel, sp::JuMP.Model, cut::Cut)
    sense = getsense(sp)

    # loop through previously visited states comparing the new cut against the
    # previous best. If it is strictly better, keep the new cut.
    push!(o.cuts, StoredCut(cut, 0))
    cut_index = length(o.cuts)
    for state in o.states
        y = cut.intercept + dot(cut.coefficients, state.state)
        if dominates(sense, y, state.best_objective)
            # if new cut is strictly better
            # decrement the counter at the old cut
            o.cuts[state.best_cut_index].non_dominated_count -= 1
            # increment the counter at the old cut
            o.cuts[cut_index].non_dominated_count += 1
            state.best_cut_index = cut_index
            state.best_objective = y
        end
    end

    # get the last state
    current_state = copy(getstage(m, ext(sp).stage).state)
    if length(current_state) == 0
        # probably in the async version
        # where we're adding a cut but haven't seen a state yet
        # or loading cuts to a new model
        return
    end

    if current_state in o.sampled_states
        return
    end
    push!(o.sampled_states, current_state)
    # now loop through the previously discovered cuts comparing them at the
    # new sampled state. If the new cut is strictly better, keep it, otherwise
    # keep the old cut
    sampled_state = SampledState(current_state,
        cut.intercept + dot(cut.coefficients, current_state),
        cut_index  # assume that the new cut is the best
    )
    push!(o.states, sampled_state)
    o.cuts[cut_index].non_dominated_count += 1

    for (i, stored_cut) in enumerate(o.cuts)
        y = stored_cut.cut.intercept + dot(stored_cut.cut.coefficients, sampled_state.state)
        if dominates(sense, y, sampled_state.best_objective)
            # if new cut is strictly better
            # decrement the counter at the old cut
            o.cuts[sampled_state.best_cut_index].non_dominated_count -= 1
            # increment the counter at the old cut
            o.cuts[i].non_dominated_count += 1
            sampled_state.best_cut_index = i
            sampled_state.best_objective = y
        end
    end
end

function validcuts(o::LevelOneCutOracle)
    active_cuts = Cut[]
    for stored_cut in o.cuts
        if stored_cut.non_dominated_count > 0
            push!(active_cuts, stored_cut.cut)
        end
    end
    active_cuts
end
