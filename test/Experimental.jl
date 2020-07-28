#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

import GLPK
import JSON
import JSONSchema
using SDDP
using Test

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
                    sum(C[i] * x[i].out for i = 1:N) -
                    S[ω] * s[ω] - s[ω] * S[ω] +
                    ω
                )
            )
        end
    end
    return model
end

download(
    "https://odow.github.io/StochOptFormat/versions/sof-0.1.schema.json",
    "sof.schema.json"
)
const SCHEMA = JSONSchema.Schema(
    JSON.parsefile("sof.schema.json"; use_mmap = false)
)

@testset "StochOptFormat" begin
    @testset "Min: Read and write to file" begin
        base_model = _create_model(true)
        set_optimizer(base_model, GLPK.Optimizer)
        SDDP.train(base_model; iteration_limit = 50, print_level = 0)

        model = _create_model(true)
        SDDP.write_to_file(model, "experimental.sof.json"; test_scenarios = 10)
        @test isvalid(JSON.parsefile("experimental.sof.json"), SCHEMA)
        set_optimizer(model, GLPK.Optimizer)
        SDDP.train(model; iteration_limit = 50, print_level = 0)

        new_model, test_scenarios = SDDP.read_from_file("experimental.sof.json")
        set_optimizer(new_model, GLPK.Optimizer)
        SDDP.train(new_model; iteration_limit = 50, print_level = 0)

        @test isapprox(
            SDDP.calculate_bound(base_model),
            SDDP.calculate_bound(model);
            atol = 1e-6
        )

        @test isapprox(
            SDDP.calculate_bound(base_model),
            SDDP.calculate_bound(new_model);
            atol = 1e-6
        )

        scenarios = SDDP.evaluate(new_model, test_scenarios)
        @test length(scenarios["problem_sha256_checksum"]) == 64
        @test length(scenarios["scenarios"]) == 10
        @test length(scenarios["scenarios"][1]) == 3
        node_1_1 = scenarios["scenarios"][1][1]
        @test isapprox(node_1_1["objective"], 9.6; atol = 1e-8)
        @test node_1_1["primal"]["d"] == 2
        @test isapprox(node_1_1["primal"]["x[1]_out"], 1; atol = 1e-8)
        @test isapprox(node_1_1["primal"]["x[2]_out"], 12; atol = 1e-8)
    end

    @testset "Max: Read and write to file" begin
        base_model = _create_model(false)
        set_optimizer(base_model, GLPK.Optimizer)
        SDDP.train(base_model; iteration_limit = 50, print_level = 0)

        model = _create_model(false)
        SDDP.write_to_file(model, "experimental.sof.json"; test_scenarios = 10)
        @test isvalid(JSON.parsefile("experimental.sof.json"), SCHEMA)
        set_optimizer(model, GLPK.Optimizer)
        SDDP.train(model; iteration_limit = 50, print_level = 0)

        new_model, test_scenarios = SDDP.read_from_file("experimental.sof.json")
        set_optimizer(new_model, GLPK.Optimizer)
        SDDP.train(new_model; iteration_limit = 50, print_level = 0)

        @test isapprox(
            SDDP.calculate_bound(base_model),
            SDDP.calculate_bound(model);
            atol = 1e-6
        )

        @test isapprox(
            SDDP.calculate_bound(base_model),
            SDDP.calculate_bound(new_model);
            atol = 1e-6
        )

        scenarios = SDDP.evaluate(new_model, test_scenarios)
        @test length(scenarios["problem_sha256_checksum"]) == 64
        @test length(scenarios["scenarios"]) == 10
        @test length(scenarios["scenarios"][1]) == 3
        node_1_1 = scenarios["scenarios"][1][1]
        @test isapprox(node_1_1["objective"], -9.6; atol = 1e-8)
        @test node_1_1["primal"]["d"] == 2
        @test isapprox(node_1_1["primal"]["x[1]_out"], 1; atol = 1e-8)
        @test isapprox(node_1_1["primal"]["x[2]_out"], 12; atol = 1e-8)
    end

    @testset "kwarg to Base.write" begin
        model = _create_model(true)
        SDDP.write_to_file(
            model,
            "experimental.sof.json";
            test_scenarios = 0,
            name = "Experimental",
            description = "Experimental model",
            author = "Oscar Dowson",
            date = "1234-56-78",
        )
        data = JSON.parsefile("experimental.sof.json", use_mmap = false)
        @test isvalid(data, SCHEMA)
        @test data["description"] == "Experimental model"
        @test data["author"] == "Oscar Dowson"
        @test data["date"] == "1234-56-78"
    end

    @testset "Error: existing cuts" begin
        model = _create_model(true)
        set_optimizer(model, GLPK.Optimizer)
        SDDP.train(model; iteration_limit = 1, print_level = 0)
        err = ErrorException(
            "StochOptFormat does not support writing after a call to " *
            "`SDDP.train`."
        )
        @test_throws err Base.write(IOBuffer(), model)
    end

    @testset "Error: belief states" begin
        model = _create_model(true; belief = true)
        err = ErrorException("StochOptFormat does not support belief states.")
        @test_throws err Base.write(IOBuffer(), model)
    end

    @testset "Error: objective states" begin
        model = _create_model(true; objective_state = true)
        err = ErrorException("StochOptFormat does not support objective states.")
        @test_throws err Base.write(IOBuffer(), model)
    end
end

if isfile("experimental.sof.json")
    rm("experimental.sof.json")
end
if isfile("sof.schema.json")
    rm("sof.schema.json")
end

@testset "slptestset" begin
    model, test_scenarios = SDDP.read_from_file(
        joinpath(dirname(@__DIR__), "examples/StochOptFormat/electric.sof.json")
    )
    set_optimizer(model, GLPK.Optimizer)
    SDDP.train(model; iteration_limit = 20, print_level = 0)
    @test isapprox(SDDP.calculate_bound(model), 381.8533; atol = 1e-3)
    scenarios = SDDP.evaluate(model, test_scenarios)
    @test length(scenarios["problem_sha256_checksum"]) == 64
    @test length(scenarios["scenarios"]) == 3
end
