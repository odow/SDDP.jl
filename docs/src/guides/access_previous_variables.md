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

## Variable lead times

What do you do when the lead times are a stochastic variable?

In this example, we consider the lead time to follow a truncated Geometric 
distribution.

### Tracking the orders
We add a dimension on the x_orders variable to track the orders in transit, assuming 
they can take 1:T lead times. Thus, x_orders[1] are arriving in the inventory and 
x_orders[2] will arrive in the next stage.

### Lead times as stochastic variable
The second trick here is to use ``\omega`` as a stochastic variable to represent the lead 
together with set_normalized_coefficient to enable u_buy on the ``\omega`` is equal to i 
and disable to the rest. For example, if ``\omega`` = 2, it means that it will take 2 stages 
for the order to arrive. In this case, we will have (for `T` =4):
```julia
c_orders[1], x_orders[2].in + 0 * u_buy == x_orders[1].out
c_orders[2], x_orders[3].in + 1 * u_buy == x_orders[2].out
c_orders[3], x_orders[4].in + 0 * u_buy == x_orders[3].out
c_orders[4], x_orders[5].in + 0 * u_buy == x_orders[4].out
```
Let us take a look at c_orders[2]. Decision variable x_orders[2].out will receive 
the orders that were already placed and will arrive in 3 stages x_orders[3].in 
plus u_buy. In the next stage, this order will pushed to x_orders[1].out and, in 
the following stage, to x_orders[1].in and into the inventory.

```@repl
using SDDP
import HiGHS
import Distributions
T = 10
model = SDDP.LinearPolicyGraph(
    stages = 20,
    sense = :Max,
    upper_bound = 1000,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variables(sp, begin
        x_inventory >= 0, SDDP.State, (initial_value = 0)
        # Add an extra dimention on the orders to track the orders lead times
        x_orders[1:T+1], SDDP.State, (initial_value = 0)
        0 <= u_buy <= 10
        u_sell >= 0
    end)
    fix(x_orders[T+1].out, 0)
    @stageobjective(sp, u_sell)
    @constraints(sp, begin
        # Shift the orders one stage 
        c_orders[i=1:T], x_orders[i+1].in + 1 * u_buy == x_orders[i].out
        # x_orders[1].in are arriving on the inventory
        x_inventory.out == x_inventory.in - u_sell + x_orders[1].in
    end)
    # 
    Ω = 1:T
    P = Distributions.pdf.(Distributions.Geometric(1 / 5), 0:T-1)
    P ./= sum(P)
    SDDP.parameterize(sp, Ω, P) do ω
        # Rewrite the constraint c_orders[i=1:T] indicating how many stages
        #  ahead the order will arrive (ω)
        #   x_orders[i+1].in + 1 * u_buy == x_orders[i].out
        # if ω == i and
        #   x_orders[i+1].in + 0 * u_buy == x_orders[i].out
        # if ω != i.
        for i in Ω
            set_normalized_coefficient(c_orders[i], u_buy, ω == i ? 1 : 0)
        end
    end
end
```
