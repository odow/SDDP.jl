#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestMSPFormat

using SDDP
using Test
import HiGHS

import SDDP: MSPFormat

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

function test_get_constant()
    @test MSPFormat._get_constant(Any[1.0]) == 1.0
    @test MSPFormat._get_constant(Any[2.4]) == 2.4
    @test MSPFormat._get_constant(Any["inf"]) == Inf
    @test MSPFormat._get_constant(Any["-inf"]) == -Inf
    state = Dict{String,Any}("inflow" => 12.0)
    @test MSPFormat._get_constant(Any[1.0], state) == 1.0
    @test MSPFormat._get_constant(Any[2.4], state) == 2.4
    @test MSPFormat._get_constant(Any["inf"], state) == Inf
    @test MSPFormat._get_constant(Any["-inf"], state) == -Inf
    terms = [Dict("ADD" => "inflow"), Dict("ADD" => 0.0)]
    @test MSPFormat._get_constant(terms) === terms
    @test MSPFormat._get_constant(terms, state) === 12.0
    terms = Any[
        Dict("ADD" => "inflow"),
        Dict("ADD" => 200.0),
        Dict("ADD" => Any[0.0]),
    ]
    @test MSPFormat._get_constant(terms) === terms
    @test MSPFormat._get_constant(terms, state) === 212.0
    terms = Any[
        Dict("ADD" => "inflow"),
        Dict("ADD" => 200.0),
        Dict("ADD" => Any[1.0]),
    ]
    @test MSPFormat._get_constant(terms) === terms
    @test MSPFormat._get_constant(terms, state) === 213.0
    return
end

function test_set_type()
    @test MSPFormat._set_type(1.0, "EQ") == JuMP.MOI.EqualTo(1.0)
    @test MSPFormat._set_type(1.2, "EQ") == JuMP.MOI.EqualTo(1.2)
    @test MSPFormat._set_type(Any[], "EQ") == JuMP.MOI.EqualTo(0.0)
    @test MSPFormat._set_type(1.0, "LEQ") == JuMP.MOI.LessThan(1.0)
    @test MSPFormat._set_type(1.2, "LEQ") == JuMP.MOI.LessThan(1.2)
    @test MSPFormat._set_type(Any[], "LEQ") == JuMP.MOI.LessThan(0.0)
    @test MSPFormat._set_type(1.0, "GEQ") == JuMP.MOI.GreaterThan(1.0)
    @test MSPFormat._set_type(1.2, "GEQ") == JuMP.MOI.GreaterThan(1.2)
    @test MSPFormat._set_type(Any[], "GEQ") == JuMP.MOI.GreaterThan(0.0)
    return
end

function test_SimpleHydroThermal()
    problem = joinpath(@__DIR__, "Simple_Hydrothermal")
    if haskey(ENV, "CI")
        @test !isfile("$problem.problem.json")
        return  # Skip tests
    end
    model = MSPFormat.read_from_file(problem)
    JuMP.set_optimizer(model, HiGHS.Optimizer)
    SDDP.train(model; iteration_limit = 10, print_level = 0)
    @test ≈(SDDP.calculate_bound(model), 8333.3333, atol = 1e-4)
    return
end

function test_SimpleHydroThermal_round_trip()
    problem = joinpath(@__DIR__, "Simple_Hydrothermal")
    if haskey(ENV, "CI")
        @test !isfile("$problem.problem.json")
        return  # Skip tests
    end
    src = MSPFormat.read_from_file(problem)
    SDDP.write_to_file(src, "$problem.sof.json")
    model, validation = SDDP.read_from_file("$problem.sof.json")
    @test validation === nothing
    JuMP.set_optimizer(model, HiGHS.Optimizer)
    SDDP.train(model; iteration_limit = 10, print_level = 0)
    @test ≈(SDDP.calculate_bound(model), 8333.3333, atol = 1e-4)
    rm("$problem.sof.json")
    return
end

end  # module

TestMSPFormat.runtests()
