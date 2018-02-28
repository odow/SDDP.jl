#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

@testset "States" begin
    rhs   = [1,2,3]
    duals = [0.5, 1.5, 2.5]
    newrhs = [1.1, 1.2, 1.3]

    @testset "Single element" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @state(m, x, x0 == rhs[1])
        @test length(subproblem.states) == 1
        @state(m, y >= 0, y0 == rhs[2])
        @test length(subproblem.states) == 2
        @state(m, 1 <= z <= 10, z0 == rhs[3])
        @test length(subproblem.states) == 3
        v = [x, y, z]
        v0 = [x0, y0, z0]
        m.linconstrDuals = duals
        for i in 1:3
            s = subproblem.states[i]
            @test s.variable == v[i]
            c = m.linconstr[s.constraint.idx]
            @test c.lb == c.ub == rhs[i]
            SDDP.setvalue!(s, newrhs[i])
            @test c.lb == c.ub == newrhs[i]
            @test c.terms == JuMP.AffExpr(v0[i])
            @test JuMP.getdual(s) == duals[i]
        end
    end

    @testset "JuMP Array" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        m.ext[:SDDP] = subproblem
        @state(m, x[i=1:3], x0 == i)
        @test length(subproblem.states) == 3
        m.linconstrDuals = duals
        for i in 1:3
            s = subproblem.states[i]
            @test s.variable == x[i]
            c = m.linconstr[s.constraint.idx]
            @test c.lb == c.ub == rhs[i]
            SDDP.setvalue!(s, newrhs[i])
            @test c.lb == c.ub == newrhs[i]
            @test c.terms == JuMP.AffExpr(x0[i])
            @test JuMP.getdual(s) == duals[i]
        end

    end

    @testset "JuMP Dict" begin
        set     = [:a, :b, :c]
        _rhs    = Dict(zip(set, rhs))
        _duals  = Dict(zip(set, duals))
        _newrhs = Dict(zip(set, newrhs))
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @state(m, x[i=set], x0 == _rhs[i])
        @test length(subproblem.states) == 3
        m.linconstrDuals = duals
        for (idx, i) in enumerate(set)
            s = subproblem.states[idx]
            @test s.variable == x[i]
            c = m.linconstr[s.constraint.idx]
            @test c.lb == c.ub == _rhs[i]
            SDDP.setvalue!(s, _newrhs[i])
            @test c.lb == c.ub == _newrhs[i]
            @test c.terms == JuMP.AffExpr(x0[i])
            @test JuMP.getdual(s) == _duals[i]
        end
    end

    @testset "@states" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @states(m, begin
            x, x0 == rhs[1]
            y >= 0, y0 == rhs[2]
            1 <= z <= 10, z0 == rhs[3]
        end)
        @test length(subproblem.states) == 3
        v = [x, y, z]
        v0 = [x0, y0, z0]
        m.linconstrDuals = duals
        for i in 1:3
            s = subproblem.states[i]
            @test s.variable == v[i]
            c = m.linconstr[s.constraint.idx]
            @test c.lb == c.ub == rhs[i]
            SDDP.setvalue!(s, newrhs[i])
            @test c.lb == c.ub == newrhs[i]
            @test c.terms == JuMP.AffExpr(v0[i])
            @test JuMP.getdual(s) == duals[i]
        end
    end
end
