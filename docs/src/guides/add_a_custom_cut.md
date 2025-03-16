# Add a custom cut

Sometimes you may want to add a set of cuts to the value function that you have
computed outside of SDDP.jl. The easiest way to achieve this is via the
[`read_cuts_from_file`](@ref) function.

To see the naming convention required by [`read_cuts_from_file`](@ref), train
your model for one iteration, and use [`write_cuts_to_file`](@ref) to write the
cuts to file:

```@repl guide_add_a_custom_cut
using SDDP, HiGHS
function create_model()
    return SDDP.LinearPolicyGraph(;
        stages = 3,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, 0 <= x <= 100, Int, SDDP.State, initial_value = 0)
        @variable(sp, 0 <= u_p <= 200, Int)
        @variable(sp, u_o >= 0, Int)
        @variable(sp, w)
        @constraint(sp, x.out == x.in + u_p + u_o - w)
        @stageobjective(sp, 100 * u_p + 300 * u_o + 50 * x.out)
        Ω = [[100.0], [100.0, 300.0], [100.0, 300.0]]
        SDDP.parameterize(ω -> JuMP.fix(w, ω), sp, Ω[t])
    end
end
model = create_model();
SDDP.train(model; iteration_limit = 1, print_level = 0)
SDDP.write_cuts_to_file(model, "cuts.json")
print(read("cuts.json", String))
```

Then create a new file containing the cut. The formula for the cut is
```
theta >= intercept + sum(coefficients[k] * (x[k] - state[k]) for k in keys(state))
```

In this example, we add the cut ``\theta \ge 55500 - 200 * (x - 10)`` to node 1:
```@jldoctest guide_add_a_custom_cut
julia> model = create_model();

julia> SDDP.calculate_bound(model)
10000.0

julia> write(
           "new_cuts.json",
           """
           [{
               "node": "1",
               "single_cuts": [{
                   "state": {"x": 10.0},
                   "coefficients": {"x": -200.0},
                   "intercept": 55500.0
               }],
               "risk_set_cuts": [],
               "multi_cuts": []
           }]
           """
       );

julia> SDDP.read_cuts_from_file(model, "new_cuts.json")

julia> SDDP.calculate_bound(model)
62500.0
```

For most purposes, you can ignore the `risk_set_cuts` and `multi_cuts` fields of
the JSON file; they are used when you are using `SDDP.MULTI_CUT` and a risk
measure.
