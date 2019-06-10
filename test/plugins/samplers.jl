#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, Test, GLPK

function test_news_vendor()
    graph = SDDP.Graph(
               :root_node,
               [:decision_node,:demand_node],
               [
                   (:root_node => :decision_node, 1.0),
                   (:decision_node => :demand_node, 1.0)
               ]);

    p = 1.0 #day ahead price
    q = 2.0 #same day price

    model = SDDP.PolicyGraph(
                graph,
                sense = :Min,
                lower_bound = 0,
                optimizer = with_optimizer(GLPK.Optimizer)) do subproblem, node

                    @variable(subproblem, stock>=0, SDDP.State, initial_value = 0)
                    @variable(subproblem, 0 <= reserve)
                    @variable(subproblem, 0 <= shortage)
                    @variable(subproblem, demand)
                    @constraint(subproblem, stock.out == stock.in + reserve + shortage - demand)

                    if node == :decision_node
                            JuMP.fix(demand,0.0)
                            JuMP.fix(shortage,0.0;force=true)
                    else
                        demand_range = collect(0:100);
                        demand_probability = ones(size(demand_range));
                        demand_probability ./= sum(demand_probability);

                        SDDP.parameterize(subproblem,demand_range,demand_probability) do d
                            JuMP.fix(demand,d)
                            JuMP.fix(reserve,0.0;force=true)
                        end
                    end

                    if node == :decision_node
                        @stageobjective(subproblem,  p*reserve);
                    else
                        @stageobjective(subproblem,  q*shortage);
                    end
                end;

    @testset "CompleteSampler" begin
        SDDP.train(model,iteration_limit=100)
        results=SDDP.simulate(model,1,[:reserve,:stock,:demand,:shortage]);
        @test mean_model(results)≈ 75 atol=5
        @test lb_model(model) ≈ 75 atol=5
    end

    # @testset "MonteCarloSampler" begin
    #      SDDP.train(model,iteration_limit=100, backward_pass_sampler=SDDP.MonteCarloSampler(100))
    #      results=SDDP.simulate(model,100,[:reserve,:stock,:demand,:shortage]);
    #      @test mean_model(results)≈ 75 atol=5
    #      @test lb_model(model) ≈ 75 atol=5
    # end
end

function mean_model(result)
    objective_values = [
               sum(stage[:stage_objective] for stage in sim) for sim in result
           ];
    return round(mean(objective_values));
end

function lb_model(model)
    return round(SDDP.calculate_bound(model))
end

test_news_vendor()
