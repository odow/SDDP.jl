#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestBackwardPassSamplingSchemes

using SDDP
using Test
using HiGHS

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_CompleteSampler()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x, ω)
        end
    end
    terms = SDDP.sample_backward_noise_terms(SDDP.CompleteSampler(), model[1])
    @test terms == model[1].noise_terms
    return
end

function test_MonteCarloSampler_1()
    model = SDDP.LinearPolicyGraph(
        stages = 1,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, [1, 3], [0.9, 0.1]) do ω
            return JuMP.set_upper_bound(x, ω)
        end
    end
    term_count = 0
    for i in 1:100
        terms = SDDP.sample_backward_noise_terms(
            SDDP.MonteCarloSampler(1),
            model[1],
        )
        @test terms[1].probability == 1.0
        if terms[1].term == model[1].noise_terms[1].term
            term_count += 1
        else
            term_count -= 1
        end
    end
    @test term_count > 20
    return
end

function test_MonteCarloSampler_100()
    model = SDDP.LinearPolicyGraph(
        stages = 1,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, [1, 3], [0.9, 0.1]) do ω
            return JuMP.set_upper_bound(x, ω)
        end
    end
    terms =
        SDDP.sample_backward_noise_terms(SDDP.MonteCarloSampler(100), model[1])
    term_count = 0
    for term in terms
        @test term.probability == 0.01
        if term.term == model[1].noise_terms[1].term
            term_count += 1
        else
            term_count -= 1
        end
    end
    @test term_count > 20
    return
end

mutable struct WithStateSampler <: SDDP.AbstractBackwardSamplingScheme
    number_of_samples::Int
end
function test_WithStateSampler()
    function sample_backward_noise_terms_with_state(
        sampler::WithStateSampler,
        node::SDDP.Node,
        state::Dict{Symbol,Float64},
    )
        if state[:x] / node.index == 1.0
            return [
                SDDP.Noise((ϵ = 3.0,), 1 / sampler.number_of_samples) for
                i in 1:sampler.number_of_samples
            ]
        elseif state[:x] / node.index == 3.0
            return [
                SDDP.Noise((ϵ = 1.0,), 1 / sampler.number_of_samples) for
                i in 1:sampler.number_of_samples
            ]
        end
    end
    model = SDDP.LinearPolicyGraph(
        stages = 5,
        lower_bound = 0.0,
        direct_mode = false,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @variable(node, ϵ)
        SDDP.parameterize(node, stage * [1, 3], [0.9, 0.1]) do ω
            return JuMP.fix(ϵ, ω)
        end
        @constraint(node, x.out == ϵ)
    end
    forward_trajectory = SDDP.forward_pass(
        model,
        SDDP.Options(model, Dict(:x => 1.0)),
        SDDP.DefaultForwardPass(),
    )
    for node_index in 1:length(forward_trajectory.scenario_path)
        state = forward_trajectory.sampled_states[node_index]
        terms = sample_backward_noise_terms_with_state(
            WithStateSampler(100),
            model[node_index],
            state,
        )
        for term in terms
            @test term.probability == 0.01
            if state[:x] / node_index == 1.0
                @test term.term.ϵ == 3.0
            elseif state[:x] / node_index == 3.0
                @test term.term.ϵ == 1.0
            end
        end
    end

    return
end

end  # module

TestBackwardPassSamplingSchemes.runtests()
