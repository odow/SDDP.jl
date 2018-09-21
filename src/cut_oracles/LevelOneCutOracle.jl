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

function store_cut(oracle::LevelOneCutOracle, model::SDDPModel,
                   subproblem::JuMP.Model, cut::Cut)
    sense = getsense(subproblem)

    # Loop through previously visited states comparing the new cut against the
    # previous best. If it is strictly better, keep the new cut.
    push!(oracle.cuts, StoredCut(cut, 0))
    cut_index = length(oracle.cuts)
    for state in oracle.states
        height = cut.intercept + dot(cut.coefficients, state.state)
        if dominates(sense, height, state.best_objective)
            # If new cut is strictly better decrement the counter at the
            # previous best.
            oracle.cuts[state.best_cut_index].non_dominated_count -= 1
            # Increment the counter at the new cut.
            oracle.cuts[cut_index].non_dominated_count += 1
            state.best_cut_index = cut_index
            state.best_objective = height
        end
    end

    # get the last state
    current_state = copy(getstage(model, ext(subproblem).stage).state)
    if length(current_state) == 0
        # This is a special case for the asynchronous algorithm where we're
        # adding a cut but haven't seen a state yet, or for the case where we're
        # loading cuts into a new model.
        return
    end

    if current_state in oracle.sampled_states
        return
    end
    push!(oracle.sampled_states, current_state)
    # Now loop through the previously discovered cuts comparing them at the new
    # sampled state. If the new cut is strictly better, keep it, otherwise keep
    # the old cut.
    sampled_state = SampledState(current_state,
        cut.intercept + dot(cut.coefficients, current_state),
        cut_index  # Assume that the new cut is the best.
    )
    push!(oracle.states, sampled_state)
    oracle.cuts[cut_index].non_dominated_count += 1

    for (index, stored_cut) in enumerate(oracle.cuts)
        height = stored_cut.cut.intercept + dot(stored_cut.cut.coefficients,
            sampled_state.state)
        if dominates(sense, height, sampled_state.best_objective)
            # If new cut is strictly better,  decrement the counter at the old
            # cut.
            oracle.cuts[sampled_state.best_cut_index].non_dominated_count -= 1
            # Increment the counter at the new cut.
            oracle.cuts[index].non_dominated_count += 1
            sampled_state.best_cut_index = index
            sampled_state.best_objective = height
        end
    end
end

function valid_cuts(oracle::LevelOneCutOracle)
    return Cut[stored_cut.cut for stored_cut in oracle.cuts
               if stored_cut.non_dominated_count > 0]
end

function all_cuts(oracle::LevelOneCutOracle)
    return [stored_cut.cut for stored_cut in oracle.cuts]
end
