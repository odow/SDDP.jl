using Kokako

model = Kokako.PolicyGraph(
            Kokako.Linear(stages = 5),
            root_node = Dict(:x => 1.0, :y => 2.0),
            optimizer = with_optimizer(GLPK.Optimizer),
            sense = :Min
                ) do subproblem, index, incoming_state, outgoing_state
    # index = 1, 2, 3, 4, 5
    @variable(subproblem, control >= 0)
    @constraints(subproblem, begin
        con, outgoing == incoming + control
        0 <= outgoing_state <= 1
    end)
    @objective(subproblem)

    # If not probability given in 3'rd argument, defaults to uniform.
    Kokako.parameterize(subproblem, [1, 2, 3], [0.5, 0.2, 0.3]) do ω
        JuMP.setRHS(subproblem, con = ω)
    end
end

markov_chain = Kokako.Markovian([0.5 0.5; 0.3 0.4], stages=5,
                                initial_probability=[1.0, 0.0])
# ... or ...
markov_chain = Kokako.Markovian([
    [1.0],
    [0.5, 0.5],
    [0.5 0.5; 0.3 0.4],
    [0.5 0.5; 0.3 0.4],
    [0.5 0.5; 0.3 0.4]
], initial_probability=[1.0])

model = Kokako.PolicyGraph(markov_chain, sense=:Min) do subproblem, index
    # index = (1,1), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (5,1), (5,2)
    @state(subproblem, 0 <= outgoing <= 1, incoming == 1)
    @variables(subproblem, begin
        control >= 0
    end)
    con = @constraint(subproblem, outgoing == incoming + control)
    @objective(subproblem)

    Kokako.parameterize(subproblem, [1, 2, 3], [0.5, 0.2, 0.3]) do ω
        JuMP.set_lower_bound(control, 1.0)
        JuMP.set_rhs(con, ω)
    end
end
