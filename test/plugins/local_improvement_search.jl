#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestLocalImprovementSearch

using Test
import SDDP: LocalImprovementSearch
import HiGHS

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

function test_x_squared()
    calls = 0
    f, x = LocalImprovementSearch.minimize([0.0]) do x
        f = 2.1 + (x[1] - 1.1)^2
        f′ = [2 * (x[1] - 1.1)]
        calls += 1
        return f, f′
    end
    @info "squared = $(calls)"
    @test isapprox(f, 2.1, atol = 1e-6)
    @test isapprox(x, [1.1], atol = 1e-4)
    return
end

function test_exp()
    calls = 0
    f, x = LocalImprovementSearch.minimize([1.0]) do x
        calls += 1
        if x[1] < 0.1 || x[1] > 20
            return nothing
        end
        return exp(x[1]), [exp(x[1])]
    end
    @info "exp = $(calls)"
    @test isapprox(f, exp(0.1), atol = 1e-2)
    @test isapprox(x, [0.1], atol = 1e-2)
    return
end

function test_piecewise()
    calls = 0
    f, x = LocalImprovementSearch.minimize([0.05]) do x
        calls += 1
        if x[1] < 0.0
            return nothing
        elseif 0.0 <= x[1] < 0.1
            return -0.1 - 1 * (x[1] - 0.0), [-1.0]
        elseif 0.1 <= x[1] < 0.4
            return -0.2 - 0.8 * (x[1] - 0.1), [-0.8]
        elseif 0.4 <= x[1] <= 1.0
            return -0.44 + 0.1 * (x[1] - 0.4), [0.1]
        else
            @assert 1.0 <= x[1]
            return nothing
        end
    end
    @info "piecewise = $(calls)"
    @test isapprox(f, -0.44, atol = 1e-3)
    @test isapprox(x, [0.4], atol = 1e-3)
    return
end

function test_x_squared_outer_approximation()
    calls = 0
    solver = LocalImprovementSearch.OuterApproximation(HiGHS.Optimizer)
    f, x = LocalImprovementSearch.minimize(solver, [0.0], 0.0) do x
        f = 2.1 + (x[1] - 1.1)^2
        f′ = [2 * (x[1] - 1.1)]
        calls += 1
        return f, f′
    end
    @info "OA(squared) = $(calls)"
    @test isapprox(f, 2.1, atol = 1e-6)
    @test isapprox(x, [1.1], atol = 1e-4)
    return
end

function test_exp_outer_approximation()
    calls = 0
    solver = LocalImprovementSearch.OuterApproximation(HiGHS.Optimizer)
    f, x = LocalImprovementSearch.minimize(solver, [1.0], 0.0) do x
        calls += 1
        if x[1] < 0.1 || x[1] > 20
            return nothing
        end
        return exp(x[1]), [exp(x[1])]
    end
    @info "OA(exp) = $(calls)"
    @test isapprox(f, exp(0.1), atol = 1e-2)
    @test isapprox(x, [0.1], atol = 1e-2)
    return
end

function test_piecewise_outer_approximation()
    calls = 0
    solver = LocalImprovementSearch.OuterApproximation(HiGHS.Optimizer)
    f, x = LocalImprovementSearch.minimize(solver, [0.05], -1.0) do x
        calls += 1
        if x[1] < 0.0
            return nothing
        elseif 0.0 <= x[1] < 0.1
            return -0.1 - 1 * (x[1] - 0.0), [-1.0]
        elseif 0.1 <= x[1] < 0.4
            return -0.2 - 0.8 * (x[1] - 0.1), [-0.8]
        elseif 0.4 <= x[1] <= 1.0
            return -0.44 + 0.1 * (x[1] - 0.4), [0.1]
        else
            @assert 1.0 <= x[1]
            return nothing
        end
    end
    @info "OA(piecewise) = $(calls)"
    @test isapprox(f, -0.44, atol = 1e-3)
    @test isapprox(x, [0.4], atol = 1e-3)
    return
end

end

TestLocalImprovementSearch.runtests()
