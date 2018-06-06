```@meta
CurrentModule = SDDP
```

## The first example

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
subproblem; and `t`, an `Int` that is a counter from `1` to the number of stages.
In this case, `t=1,2,3`. The keyword arguments to `SDDPModel` should be obvious.

Now we need to define build each subproblem using a mix of JuMP and SDDP.jl
syntax.

### State variables

There is one state variable in our model: `v`, the quantity of water in the
reservoir at the end of stage `t`.

```julia
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 50)
```

### Control variables

We define some JuMP variables:
```julia
@variables(sp, begin
    0 <= thermal_generation[1:2] <= 100
    0 <= hydro_generation        <= 150
         hydro_spill             >= 0
 end)
```

### Constraints


```julia
inflow = [20.0, 20.0, 10.0]
```

First, we have the water balance constraint: the volume of water at
the end of the stage must equal the volume of water at the start of the stage,
plus any inflows, less that used for generation or spilled down the spillway.
```julia
@constraint(sp,
    incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume
)
```

There is also a constraint that total generation must equal demand:
```julia
@constraint(sp,
    sum(thermal_generation) + hydro_generation == 150
)
```

### The stage objective    

Finally, there is a cost on thermal generation. The first 100 MWh are charged at
$100/MWh and the second 100 MWh are charged at $1000/MWh:
```julia
@stageobjective(sp, 100*thermal_generation[1] + 1000*thermal_generation[2] )
```
