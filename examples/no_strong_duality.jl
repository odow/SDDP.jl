using SDDP, GLPK

graph = SDDP.Graph(
    :root,
    [:node],
    [(:root => :node, 1.0), (:node => :node, 0.5)]
)

model = SDDP.PolicyGraph(
        graph, optimizer = with_optimizer(GLPK.Optimizer), 
        lower_bound = 0.0) do sp, t
    @variable(sp, x, SDDP.State, initial_value = 1.0)
    @stageobjective(sp, x.out)
    @constraint(sp, x.in == x.out)
end

SDDP.train(model, iteration_limit = 20)