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

```julia
SDDP.LinearPolicyGraph(
    stages = 10,
    sense = :Max,
    lower_bound = -100.0,
) do sp, t
    # Capacity of the generator. Decided in the first stage.
    @variable(sp, capcity, SDDP.State, initial_value = 0)
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

This is often useful if have some inventory problem with a lead-time on orders.

```julia
SDDP.LinearPolicyGraph(
    stages = 10,
    sense = :Max,
    upper_bound = 100,
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
