# Tutorial One: first steps

Hydrothermal scheduling is the most common application of stochastic dual
dynamic programming. To illustrate some of the basic functionality of SDDP.jl,
we implement a very simple model of the hydrothermal scheduling problem.

In this model, there are two generators: a thermal generator, and a hydro
generator. The thermal generator has a short-run marginal cost of \\\$50/MWh in
the first stage, \\\$100/MWh in the second stage, and \\\$150/MWh in the third
stage. The hydro generator has a short-run marginal cost of \\\$0/MWh.

We consider the problem of scheduling the generation over three time periods in
order to meet a known demand of 150 MWh in each period.

The hydro generator draws water from a reservoir which has a maximum capacity of
200 units. We assume that at the start of the first time period, the reservoir
is full. In addition to the ability to generate electricity by passing water
through the hydroelectric turbine, the hydro generator can also spill water down
a spillway (bypassing the turbine) in order to prevent the water from
over-topping the dam. We assume that there is no cost of spillage.

The objective of the optimization is to minimize the expected cost of generation
over the three time periods.

## Formulating the problem

First, we need to load some packages. For this example, we are going to use the
[Clp.jl](https://github.com/JuliaOpt/Clp.jl) package; however, you are free to
use any solver that you could normally use with JuMP.
```julia
using SDDP, JuMP, Clp
```

Next, we need to initialize our model. In our example, we are minimizing, there
are three stages, and we know a lower bound of `0.0`. Therefore, we can
initialize our model using the [`SDDPModel`](@ref) constructor:
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
example, the least-cost solution is to meet demand entirely from the hydro
generator, incurring a cost of \\\$0/MWh. Therefore, we set
`objective_bound=0.0`.

Now we need to define build each subproblem using a mix of JuMP and SDDP.jl
syntax.

### State variables

There is one state variable in our model: the quantity of water in the reservoir
at the end of stage `t`. Two add this state variable to the model, SDDP.jl
defines the [`@state`](@ref) macro.  This macro takes three arguments:
1. `sp` - the JuMP model;
2. an expression for the outgoing state variable; and
3. an expression for the incoming state variable.

The 2ⁿᵈ argument can be any valid JuMP `@variable` syntax and can include, for
example, upper and lower bounds. The 3ʳᵈ argument must be the name of the
incoming state variable, followed by `==`, and then the value of the state
variable at the root node of the policy graph. For our hydrothermal example, the
state variable can be constructed as:
```julia
@state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
```

### Control variables

We now need to define some control variables. In SDDP.jl, control variables are
just normal JuMP variables. Therefore, we can define the three variables in the
hydrothermal scheduling problem (thermal generation, hydro generation, and the
quantity of water to spill) as follows:
```julia
@variables(sp, begin
    thermal_generation >= 0
    hydro_generation   >= 0
    hydro_spill        >= 0
 end)
```

### Constraints

Before we specify the constraints, we need to create some data. For this
problem, we need the inflow to the reservoir in each stage `t=1, 2, 3`.
Therefore, we create the vector:
```julia
inflow = [50.0, 50.0, 50.0]
```
The inflow in stage `t` can be accessed as `inflow[t]`.

First, we have the water balance constraint: the volume of water at
the end of the stage must equal the volume of water at the start of the stage,
plus any inflows, less that used for generation or spilled down the spillway.
```julia
@constraint(sp,
    incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume
)
```
Note that we use `t` defined by the `SDDPModel` constructor. There is also a
constraint that total generation must equal demand of 150 MWh:
```julia
@constraint(sp,
    thermal_generation + hydro_generation == 150
)
```

### The stage objective    

Finally, there is a cost on thermal generation of \\\$50/MWh in the first stage,
\\\$100/MWh in the second stage, and \\\$150/MWh in the third stage. To add
the stage-objective, we use the aptly named `@stageobjective` macro provided by
SDDP.jl:
```julia
if t == 1
    @stageobjective(sp,  50.0 * thermal_generation )
elseif t == 2
    @stageobjective(sp, 100.0 * thermal_generation )
elseif t == 3
    @stageobjective(sp, 150.0 * thermal_generation )
end
```
!!! info
    `if` statements can be used more broadly in the subproblem definition to
    conditionally and variables and constraints into different subproblems.

We can also implement the stage-objective more succinctly using a vector:
```julia
fuel_cost = [50.0, 100.0, 150.0]
@stageobjective(sp, fuel_cost[t] * thermal_generation )
```

## Solving the problem

Putting all that we have discussed above together, we get:
```julia
using SDDP, JuMP, Clp
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
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end
```

To solve this problem, we use the [`solve`](@ref) method:

```julia
status = solve(m; max_iterations=5)
```

The return value `status` is a symbol describing why the SDDP algorithm
terminated. In this case, the value is `:max_iterations`. We discuss other
arguments to the [`solve`](@ref) method and other possible values for `status`
in future sections of this manual.

During the solve, the following log is printed to the screen.
```
-------------------------------------------------------------------------------
                          SDDP.jl © Oscar Dowson, 2017-2018
-------------------------------------------------------------------------------
    Solver:
        Serial solver
    Model:
        Stages:         3
        States:         1
        Subproblems:    3
        Value Function: Default
-------------------------------------------------------------------------------
              Objective              |  Cut  Passes    Simulations   Total
     Simulation       Bound   % Gap  |   #     Time     #    Time    Time
-------------------------------------------------------------------------------
       15.000K         5.000K        |     1    0.0      0    0.0    0.0
        5.000K         5.000K        |     2    0.0      0    0.0    0.0
        5.000K         5.000K        |     3    0.0      0    0.0    0.0
        5.000K         5.000K        |     4    0.0      0    0.0    0.0
        5.000K         5.000K        |     5    0.0      0    0.0    0.0
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         5
        Termination Status: max_iterations
===============================================================================
```

The header and footer of the output log contain self-explanatory statistics
about the problem. The numeric columns are worthy of description. Each row
corresponds to one iteration of the SDDP algorithm.

The left half of the log relates to the objective of the problem. In the
*Simulation* column, we give the cumulative cost of each forward pass. In the
*Bound* column, we give the lower bound (upper if maximizing) obtained after the
backward pass has completed in each iteration. Ignore the *% Gap* column for
now, that is addressed in [Tutorial Two: RHS noise](@ref).

The right half of the log displays timing statistics. *Cut Passes* displays the
number of cutting iterations conducted (in *#*) and the time it took to (in
*Time*). Ignore the *Simulations* columns for now, they are addressed in
Tutorial [Tutorial Two: RHS noise](@ref). Finally, the *Total Time* column
records the total time spent solving the problem.

This log can be silenced by setting the `print_level` keyword argument to
`solve` to `0`. In addition, the log will be written to the file given by the
`log_file` keyword argument (this is off by default).

## Understanding the solution

The first thing we want to do is to query the lower (upper if maximizing) bound
of the solution. This can be done via the [`getbound`](@ref) function:
```julia
getbound(m)
```
This returns the value of the *Bound* column in the last row in the output table
above. In this example, the bound is `5000.0`.

Then, we can perform a Monte Carlo simulation of the policy using the
[`simulate`](@ref) function. It takes three arguments. The first is the
[`SDDPModel`](@ref) `m`. The second is the number of replications to perform.
The third is a vector of variable names to record the value of at each stage and
replication. Since our example is deterministic, it is sufficient to perform a
single replication:
```julia
simulation_result = simulate(m,
    1,
    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill]
)
```
The return value, `simulation_result`, is a vector of dictionaries containing
one element for each Monte Carlo replication. In this case,
`length(simulation_result) = 1`. The keys of the dictionary are the variable
symbols given in the `simulate` function, and their associated values are
vectors, with one element for each stage, or the variable value in the simulated
solution. For example, we can query the optimal quantity of hydro generation in
each stage as follows:
```julia
julia> simulation_result[1][:hydro_generation]
3-element Array{Any, 1}:
  50.0
 150.0
 150.0
```

This concludes our first very simple tutorial for SDDP.jl. In the next tutorial,
[Tutorial Two: RHS noise](@ref), we introduce stagewise-independent noise into
the model.
