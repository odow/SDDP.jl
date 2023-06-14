using SDDP, Gurobi

function main(capacity_cost)
    graph = SDDP.LinearGraph(2)
    SDDP.add_edge(graph, 2 => 2, 0.95)
    model = SDDP.PolicyGraph(
        graph;
        sense = :Min,
        lower_bound = 0.0,
        optimizer = Gurobi.Optimizer,
    ) do sp, node
        @variable(sp, 0 <= x_capacity <= 400, SDDP.State, initial_value = 0)
        @variable(sp, 0 <= x_prod, SDDP.State, initial_value = 0)
        if node == 1
            @stageobjective(sp, capacity_cost * x_capacity.out)
            @constraint(sp, x_prod.out == x_prod.in)
        else
            @variable(sp, 0 <= u_prod <= 200)
            @variable(sp, u_overtime >= 0)
            @stageobjective(sp, 100u_prod + 300u_overtime + 50x_prod.out)
            @constraint(sp, x_capacity.out == x_capacity.in)
            @constraint(sp, x_prod.out <= x_capacity.in)
            @constraint(sp, c_bal, x_prod.out == x_prod.in + u_prod + u_overtime)
            SDDP.parameterize(sp,  [100, 300]) do ω
                JuMP.set_normalized_rhs(c_bal, -ω)
            end
        end
        return
    end
    SDDP.train(
        model;
        forward_pass = SDDP.RegularizedForwardPass(),
    )
    results = SDDP.simulate(model, 1, [:x_capacity])
    return results[1][1][:x_capacity].out
end

# results = main.([0, 100, 200, 300, 400])

main(100)
