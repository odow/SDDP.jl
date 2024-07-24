#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestExperimental

using SDDP
using Test
import Downloads
import HiGHS
import JSON
import JSONSchema

Downloads.download(
    "https://odow.github.io/StochOptFormat/schemas/sof-1.schema.json",
    "sof.schema.json",
)
const SCHEMA =
    JSONSchema.Schema(JSON.parsefile("sof.schema.json"; use_mmap = false))

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    if isfile("experimental.sof.json")
        rm("experimental.sof.json")
    end
    if isfile("sof.schema.json")
        rm("sof.schema.json")
    end
    return
end

function _create_model(
    minimization::Bool;
    belief::Bool = false,
    objective_state::Bool = false,
)
    graph = SDDP.LinearGraph(3)
    if belief
        SDDP.add_ambiguity_set(graph, [1])
        SDDP.add_ambiguity_set(graph, [2, 3])
    end
    model = SDDP.PolicyGraph(
        graph,
        sense = minimization ? :Min : :Max,
        lower_bound = -50.0,
        upper_bound = 50.0,
    ) do sp, t
        N = 2
        C = [0.2, 0.7]
        S = 2 .+ [0.33, 0.54]
        DEMAND = [2, 10]
        if objective_state
            SDDP.add_objective_state(
                (y, w) -> y + ω,
                sp,
                initial_value = 0.0,
                lower_bound = 0.0,
                upper_bound = 1.0,
                lipschitz = 1.0,
            )
        end
        @variable(sp, x[1:N] >= 0, SDDP.State, initial_value = 0.0)
        @variables(sp, begin
            s[i = 1:N] >= 0
            d
        end)
        @constraints(sp, begin
            [i = 1:N], s[i] <= x[i].in
            c, sum(s) <= d + 1
        end)
        SDDP.parameterize(sp, t == 1 ? [1] : 1:length(DEMAND)) do ω
            JuMP.fix(d, DEMAND[ω])
            set_upper_bound(s[1], 0.1 * ω)
            set_lower_bound(x[1].out, ω)
            set_normalized_rhs(c, ω)
            sgn = minimization ? 1.0 : -1.0
            @stageobjective(
                sp,
                sgn * (
                    sum(C[i] * x[i].out for i in 1:N) - S[ω] * s[ω] -
                    s[ω] * S[ω] + ω
                )
            )
        end
    end
    return model
end

function test_write_to_file_Min()
    base_model = _create_model(true)
    set_optimizer(base_model, HiGHS.Optimizer)
    SDDP.train(base_model; iteration_limit = 50, print_level = 0)
    model = _create_model(true)
    SDDP.write_to_file(
        model,
        "experimental.sof.json";
        validation_scenarios = 10,
        sampling_scheme = SDDP.PSRSamplingScheme(2),
    )
    set_optimizer(model, HiGHS.Optimizer)
    SDDP.train(model; iteration_limit = 50, print_level = 0)
    new_model, validation_scenarios =
        SDDP.read_from_file("experimental.sof.json")
    set_optimizer(new_model, HiGHS.Optimizer)
    SDDP.train(new_model; iteration_limit = 50, print_level = 0)
    @test isapprox(
        SDDP.calculate_bound(base_model),
        SDDP.calculate_bound(model);
        atol = 1e-6,
    )
    @test isapprox(
        SDDP.calculate_bound(base_model),
        SDDP.calculate_bound(new_model);
        atol = 1e-6,
    )
    scenarios = SDDP.evaluate(new_model, validation_scenarios)
    @test length(scenarios["problem_sha256_checksum"]) == 64
    @test length(scenarios["scenarios"]) == 10
    @test length(scenarios["scenarios"][1]) == 3
    node_1_1 = scenarios["scenarios"][1][1]
    @test isapprox(node_1_1["objective"], 9.6; atol = 1e-8)
    @test node_1_1["primal"]["d"] == 2
    @test isapprox(node_1_1["primal"]["x[1]_out"], 1; atol = 1e-8)
    @test isapprox(node_1_1["primal"]["x[2]_out"], 12; atol = 1e-8)
    demands = map(scenarios["scenarios"]) do s
        return [si["primal"]["d"] for si in s]
    end
    for i in 3:2:10
        @test demands[1] == demands[i]
        @test demands[2] == demands[i+1]
    end
    return
end

function test_write_to_file_Max()
    base_model = _create_model(false)
    set_optimizer(base_model, HiGHS.Optimizer)
    SDDP.train(base_model; iteration_limit = 50, print_level = 0)
    model = _create_model(false)
    SDDP.write_to_file(
        model,
        "experimental.sof.json";
        validation_scenarios = 10,
    )
    set_optimizer(model, HiGHS.Optimizer)
    SDDP.train(model; iteration_limit = 50, print_level = 0)
    new_model, validation_scenarios =
        SDDP.read_from_file("experimental.sof.json")
    set_optimizer(new_model, HiGHS.Optimizer)
    SDDP.train(new_model; iteration_limit = 50, print_level = 0)
    @test isapprox(
        SDDP.calculate_bound(base_model),
        SDDP.calculate_bound(model);
        atol = 1e-6,
    )
    @test isapprox(
        SDDP.calculate_bound(base_model),
        SDDP.calculate_bound(new_model);
        atol = 1e-6,
    )
    scenarios = SDDP.evaluate(new_model, validation_scenarios)
    @test length(scenarios["problem_sha256_checksum"]) == 64
    @test length(scenarios["scenarios"]) == 10
    @test length(scenarios["scenarios"][1]) == 3
    node_1_1 = scenarios["scenarios"][1][1]
    @test isapprox(node_1_1["objective"], -9.6; atol = 1e-8)
    @test node_1_1["primal"]["d"] == 2
    @test isapprox(node_1_1["primal"]["x[1]_out"], 1; atol = 1e-8)
    @test isapprox(node_1_1["primal"]["x[2]_out"], 12; atol = 1e-8)
end

function test_write_kwarg()
    model = _create_model(true)
    SDDP.write_to_file(
        model,
        "experimental.sof.json";
        validation_scenarios = 0,
        name = "Experimental",
        description = "Experimental model",
        author = "Oscar Dowson",
        date = "1234-56-78",
    )
    data = JSON.parsefile("experimental.sof.json", use_mmap = false)
    @test data["description"] == "Experimental model"
    @test data["author"] == "Oscar Dowson"
    @test data["date"] == "1234-56-78"
    return
end

function test_error_existing_cuts()
    model = _create_model(true)
    set_optimizer(model, HiGHS.Optimizer)
    SDDP.train(model; iteration_limit = 1, print_level = 0)
    err = ErrorException(
        "StochOptFormat does not support writing after a call to " *
        "`SDDP.train`.",
    )
    @test_throws err Base.write(IOBuffer(), model)
    return
end

function test_error_belief_states()
    model = _create_model(true; belief = true)
    err = ErrorException("StochOptFormat does not support belief states.")
    @test_throws err Base.write(IOBuffer(), model)
    return
end

function test_error_objective_states()
    model = _create_model(true; objective_state = true)
    err = ErrorException("StochOptFormat does not support objective states.")
    @test_throws err Base.write(IOBuffer(), model)
    return
end

function test_slptestset()
    model, validation_scenarios =
        SDDP.read_from_file(joinpath(@__DIR__, "electric.sof.json"))
    set_optimizer(model, HiGHS.Optimizer)
    SDDP.train(model; iteration_limit = 20, print_level = 0)
    @test isapprox(SDDP.calculate_bound(model), 381.8533; atol = 1e-3)
    scenarios = SDDP.evaluate(model, validation_scenarios)
    @test length(scenarios["problem_sha256_checksum"]) == 64
    @test length(scenarios["scenarios"]) == 3
    return
end

end  # module

TestExperimental.runtests()
