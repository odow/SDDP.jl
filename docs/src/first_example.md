```@meta
CurrentModule = SDDP
```

## Example: simple hydrothermal scheduling

Hydrothermal scheduling is the most common application of stochastic dual
dynamic programming. To illustrate some of the basic functionality of SDDP.jl,
we implement a very simple model of the hydrothermal scheduling problem.

In this model, where are two generators: a thermal generator (with a cost of
$100/MWh), and a hydro-generator (with a cost of $0/MWh). We consider the
problem of scheduling the generation over three time periods in order to meet a
known demand of 150 MWh in each period.

The hydro-generator draws water from a reservoir which has a maximum capacity of
200 units. We assume that at the start of the first time period, the reservoir
is full. In addition to the ability to generate electricity by passing water
through the hydroelectric turbine, the hydro-generator can also spill water down
a spillway (bypassing the turbine) in order to prevent the water from
over-topping the dam. We assume that there is no cost of spillage.

The objective of the optimization is to minimize the expected cost of generation
over the three time periods.

### Initialization

First, we need to load some packages. For this example, we are going to use the
[Clp.jl](https://github.com/JuliaOpt/Clp.jl) package; however, you are free to
use any solver that you could normally use with JuMP.
```julia
using SDDP, JuMP, Clp
```

Next, we need to initialize our model. In our example, we are minimizing, there
are three stages, and we know a lower bound of `0.0`. Therefore, we can
initialize our model using the `SDDPModel` constructor:
```julia
m = SDDPModel(
                  sense = :Min,          
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0
                                        ) do sp, t
    # ... stuff to go here ...
end
```
If you haven't seen the `do sp, t ... end` syntax before, this syntax is
equivalent to the following:
```julia
function define_subproblem(sp::JuMP.Model, t::Int)
    # ... stuff to go here ...
end
m = SDDPModel(
    define_subproblem,
    sense           = :Min,          
    stages          = 3,
    solver          = ClpSolver(),
    objective_bound = 0.0
)
```
The function `define_subproblem` (although you can call it anything you like)
takes two arguments: `sp`, a `JuMP.Model` that we will use to build each
subproblem; and `t`, an `Int` that is a counter from `1` to the number of
stages. In this case, `t=1, 2, 3`. The `sense`, `stages`, and `solver` keyword
arguments to `SDDPModel` should be obvious; however, the `objective_bound` is
worth explaining.

In order to solve a model using SDDP, we need to define a valid lower bound for
every subproblem. (See [Introduction to SDDP](@ref) for details.) In this
example, the least-cost solution is to meet demand entirely from the
hydro-generator, incurring a cost of $0. Therefore, we set
`objective_bound=0.0`.

Now we need to define build each subproblem using a mix of JuMP and SDDP.jl
syntax.

### State variables

There is one state variable in our model: the quantity of water in the
reservoir at the end of stage `t`.

```julia
@state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
```

### Control variables

We define some JuMP variables:
```julia
@variables(sp, begin
    thermal_generation >= 0
    hydro_generation   >= 0
    hydro_spill        >= 0
 end)
```

### Constraints


```julia
inflow = [50.0, 50.0, 50.0]
```

First, we have the water balance constraint: the volume of water at
the end of the stage must equal the volume of water at the start of the stage,
plus any inflows, less that used for generation or spilled down the spillway.
```julia
@constraint(sp,
    incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume
)
```

There is also a constraint that total generation must equal demand of 150 MWh:
```julia
@constraint(sp,
    thermal_generation + hydro_generation == 150
)
```

### The stage objective    

Finally, there is a cost on thermal generation of $100/MWh:
```julia
@stageobjective(sp, 100 * thermal_generation )
```

### All together

Putting all that we have discussed above together, we
get:
```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0
                                        ) do sp, t
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
     end)
    inflow = [50.0, 50.0, 50.0]
    @constraints(sp, begin
        incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume
        thermal_generation + hydro_generation == 150
    end)
    @stageobjective(sp, 100 * thermal_generation )
end
```
