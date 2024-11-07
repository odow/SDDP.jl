# Access variables from a previous stage

A common question is "how do I use a variable from a previous stage in a
constraint?"

!!! info
    If you want to use a variable from a previous stage, it must be a state
    variable.

Here are some examples:

## Access a first-stage decision in a future stage

This is often useful if your first-stage decisions are capacity-expansion type
decisions (e.g., you choose first how much capacity to add, but because it takes
time to build, it only shows up in some future stage).

```@repl
using SDDP, HiGHS
SDDP.LinearPolicyGraph(
    stages = 10,
    sense = :Max,
    upper_bound = 100.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    # Capacity of the generator. Decided in the first stage.
    @variable(sp, capacity >= 0, SDDP.State, initial_value = 0)
    # Quantity of water stored.
    @variable(sp, reservoir >= 0, SDDP.State, initial_value = 0)
    # Quantity of water to use for electricity generation in current stage.
    @variable(sp, generation >= 0)
    if t == 1
        # There are no constraints in the first stage, but we need to push the
        # initial value of the reservoir to the next stage.
        @constraint(sp, reservoir.out == reservoir.in)
        # Since we're maximizing profit, subtract cost of capacity.
        @stageobjective(sp, -capacity.out)
    else
        # Water balance constraint.
        @constraint(sp, balance, reservoir.out - reservoir.in + generation == 0)
        # Generation limit.
        @constraint(sp, generation <= capacity.in)
        # Push capacity to the next stage.
        @constraint(sp, capacity.out == capacity.in)
        # Maximize generation.
        @stageobjective(sp, generation)
        # Random inflow in balance constraint.
        SDDP.parameterize(sp, rand(4)) do w
            set_normalized_rhs(balance, w)
        end
    end
end
```

## Access a decision from N stages ago

This is often useful if you have some inventory problem with a lead time on orders.
In the code below, we assume that the product has a lead time of 5 stages, and we 
use a state variable to track the decisions on the production for the last 5 stages. 
The decisions are passed to the next stage by shifting them by one stage.

```@repl
using SDDP, HiGHS
SDDP.LinearPolicyGraph(
    stages = 10,
    sense = :Max,
    upper_bound = 100,
    optimizer = HiGHS.Optimizer,
) do sp, t
    # Current inventory on hand.
    @variable(sp, inventory >= 0, SDDP.State, initial_value = 0)
    # Inventory pipeline.
    #   pipeline[1].out are orders placed today.
    #   pipeline[5].in are orders that arrive today and can be added to the
    #     current inventory.
    #   Stock moves up one slot in the pipeline each stage.
    @variable(sp, pipeline[1:5], SDDP.State, initial_value = 0)
    # The number of units to order today.
    @variable(sp, 0 <= buy <= 10)
    # The number of units to sell today.
    @variable(sp, sell >= 0)
    # Buy orders get placed in the pipeline.
    @constraint(sp, pipeline[1].out == buy)
    # Stock moves up one slot in the pipeline each stage.
    @constraint(sp, [i=2:5], pipeline[i].out == pipeline[i-1].in)
    # Stock balance constraint.
    @constraint(sp, inventory.out == inventory.in - sell + pipeline[5].in)
    # Maximize quantity of sold items.
    @stageobjective(sp, sell)
end
```

!!! warning
    You must initialize the same number of state variables in every stage, even
    if they are not used in that stage.

## Stochastic lead times

Stochastic lead times can be modeled by adding stochasticity to the pipeline
balance constraint.

The trick is to use the random variable ``\omega`` to represent the lead time,
together with `JuMP.set_normalized_coefficient` to add `u_buy` to the `i`
pipeline balance constraint when ``\omega`` is equal to `i`. For example, if
``\omega = 2`` and `T = 4`, we would have constraints:
```julia
c_pipeline[1], x_pipeline[1].out == x_pipeline[2].in + 0 * u_buy
c_pipeline[2], x_pipeline[2].out == x_pipeline[3].in + 1 * u_buy
c_pipeline[3], x_pipeline[3].out == x_pipeline[4].in + 0 * u_buy
c_pipeline[4], x_pipeline[4].out == x_pipeline[5].in + 0 * u_buy
```

```@repl
using SDDP
import HiGHS
T = 10
model = SDDP.LinearPolicyGraph(
    stages = 20,
    sense = :Max,
    upper_bound = 1000,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variables(sp, begin
        x_inventory >= 0, SDDP.State, (initial_value = 0)
        x_pipeline[1:T+1], SDDP.State, (initial_value = 0)
        0 <= u_buy <= 10
        u_sell >= 0
    end)
    fix(x_pipeline[T+1].out, 0)
    @stageobjective(sp, u_sell)
    @constraints(sp, begin
        # Shift the orders one stage 
        c_pipeline[i=1:T], x_pipeline[i].out == x_pipeline[i+1].in + 1 * u_buy
        # x_pipeline[1].in are arriving on the inventory
        x_inventory.out == x_inventory.in - u_sell + x_pipeline[1].in
    end)
    SDDP.parameterize(sp, 1:T) do ω
        # Rewrite the constraint c_pipeline[i=1:T] indicating how many stages
        # ahead the order will arrive (ω)
        # if ω == i:
        #   x_pipeline[i+1].in + 1 * u_buy == x_pipeline[i].out
        # else:
        #   x_pipeline[i+1].in + 0 * u_buy == x_pipeline[i].out
        for i in 1:T
            set_normalized_coefficient(c_pipeline[i], u_buy, ω == i ? 1 : 0)
        end
    end
end
```
