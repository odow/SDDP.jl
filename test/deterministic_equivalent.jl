#  Copyright (c) 2017-22, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestDeterministicEquivalent

using SDDP
using Test
import GLPK

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

function test_acyclic_linear()
    graph = SDDP.LinearGraph(2)
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(sp, x.out)
    end
    @test SDDP.is_cyclic(model) == false
    @test typeof(SDDP.deterministic_equivalent(model)) == JuMP.Model
    return
end

function test_cyclic_linear()
    graph = SDDP.LinearGraph(2)
    SDDP.add_edge(graph, 2 => 1, 0.9)
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(sp, x.out)
    end
    @test SDDP.is_cyclic(model) == true
    @test_throws(
        ErrorException(
            "Unable to formulate deterministic equivalent: " *
            "Cyclic policy graph detected!",
        ),
        SDDP.deterministic_equivalent(model)
    )
    return
end

function test_cyclic_single_node()
    graph = SDDP.Graph(
        :root,
        [:node],
        [(:root => :node, 1.0), (:node => :node, 0.9)],
    )
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(sp, x.out)
    end
    @test SDDP.is_cyclic(model) == true
    @test_throws(
        ErrorException(
            "Unable to formulate deterministic equivalent: " *
            "Cyclic policy graph detected!",
        ),
        SDDP.deterministic_equivalent(model)
    )
    return
end

function test_acyclic_Markovian()
    model = SDDP.MarkovianPolicyGraph(
        transition_matrices = [[0.5 0.5], [0.2 0.8; 0.8 0.2]],
        lower_bound = 0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(sp, x.out)
    end
    @test SDDP.is_cyclic(model) == false
    @test typeof(SDDP.deterministic_equivalent(model)) == JuMP.Model
    return
end

function test_cyclic_Markovian()
    graph = SDDP.MarkovianGraph([[0.5 0.5], [0.2 0.8; 0.8 0.2]])
    SDDP.add_edge(graph, (2, 1) => (1, 1), 0.9)
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(sp, x.out)
    end
    @test SDDP.is_cyclic(model) == true
    @test_throws(
        ErrorException(
            "Unable to formulate deterministic equivalent: " *
            "Cyclic policy graph detected!",
        ),
        SDDP.deterministic_equivalent(model)
    )
    return
end

function test_time_limit()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(sp, x.out)
    end
    @test_throws(
        ErrorException(
            "Unable to formulate deterministic equivalent: Time limit exceeded!",
        ),
        # We use a negative time limit to force error.
        SDDP.deterministic_equivalent(model; time_limit = -10.0)
    )
    return
end

function test_objective_states()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        SDDP.add_objective_state(
            sp,
            initial_value = 0.0,
            lower_bound = 0.0,
            upper_bound = 10.0,
            lipschitz = 10.0,
        ) do p, ω
            return p + ω
        end
        SDDP.parameterize(sp, [1, 2]) do ω
            p = SDDP.objective_state(sp)
            @stageobjective(sp, p * x.out)
        end
    end
    @test_throws(
        ErrorException(
            "Unable to formulate deterministic equivalent: Objective states detected!",
        ),
        SDDP.deterministic_equivalent(model)
    )
    return
end

function test_belief_states()
    graph = SDDP.MarkovianGraph([[0.5 0.5], [0.2 0.8; 0.8 0.2]])
    SDDP.add_ambiguity_set(graph, [(1, 1), (1, 2)])
    SDDP.add_ambiguity_set(graph, [(2, 1), (2, 2)])
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(sp, x.out)
    end
    @test_throws(
        ErrorException(
            "Unable to formulate deterministic equivalent: Belief states detected!",
        ),
        SDDP.deterministic_equivalent(model)
    )
end

function test_existing_policy()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(sp, x.out)
    end
    SDDP.train(model; iteration_limit = 2, print_level = 0)
    @test_throws(
        ErrorException(
            "Unable to formulate deterministic equivalent: Model has been used " *
            "for training. Can only form deterministic equivalent on a fresh model.",
        ),
        SDDP.deterministic_equivalent(model)
    )
    return
end

function test_constant_objective()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(sp, 1.0)
    end
    d = SDDP.deterministic_equivalent(model, GLPK.Optimizer)
    optimize!(d)
    @test objective_value(d) == 2.0
    return
end

function test_constraint_with_no_terms()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @constraint(sp, x.out <= x.out)
        @stageobjective(sp, 1.0)
    end
    d = SDDP.deterministic_equivalent(model, GLPK.Optimizer)
    optimize!(d)
    @test objective_value(d) == 2.0
    return
end

function test_quadratic_expr()
    model = SDDP.LinearPolicyGraph(stages = 2, lower_bound = 0.0) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @constraint(sp, x.in^2 <= x.out)
        @stageobjective(sp, x.out)
    end
    d = SDDP.deterministic_equivalent(model)
    @test in(
        (GenericQuadExpr{Float64,VariableRef}, MOI.LessThan{Float64}),
        list_of_constraint_types(d),
    )
    return
end

function test_quadratic_expr_no_quad_terms()
    model = SDDP.LinearPolicyGraph(stages = 2, lower_bound = 0.0) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @constraint(sp, x.in^2 <= x.out + x.in^2)
        @stageobjective(sp, x.out)
    end
    d = SDDP.deterministic_equivalent(model)
    @test in(
        (GenericQuadExpr{Float64,VariableRef}, MOI.LessThan{Float64}),
        list_of_constraint_types(d),
    )
    return
end

function test_vector_valued_functions()
    model = SDDP.LinearPolicyGraph(stages = 2, lower_bound = 0.0) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @constraint(sp, [x.in, x.out] in MOI.SOS1([1.0, 2.0]))
        @stageobjective(sp, x.out)
    end
    d = SDDP.deterministic_equivalent(model)
    @test (Vector{VariableRef}, MOI.SOS1{Float64}) in
          list_of_constraint_types(d)
    return
end

end  # module

TestDeterministicEquivalent.runtests()
