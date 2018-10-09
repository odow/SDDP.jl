# Tutorial Eleven: infinite-horizon SDDP
This tutorial discusses the use of infinite-horizon stochastic dynamic dual programming (infinite-horizon SDDP). Infinite-horizon SDDP is a methodology to find the optimal steady state policy of a multi-stage stochastic problem.

We implemented infinite-horizon SDDP using the average-cost method. Another possible method could be using the discounted-cost method however this method converges slower (than the average-cost method)

My Honors [thesis](https://github.com/shasafoster/SDDP.jl/blob/master/docs/src/assets/foster_thesis.pdf) may be useful for further understanding of the underlying theory of infinite-horizon SDDP. Ben Fulton applied infinite-horizon SDDP when modelling various scenarios in the New Zealand electricity market thus his [thesis](https://github.com/shasafoster/SDDP.jl/blob/master/docs/src/assets/fulton_thesis.pdf) may also be of interest. 

The differences between infinite-horizon SDDP and standard SDDP will be explained by formulating and solving a multistage stochastic dynamic program with each method. 

## The Problem

Recall that our model for the hydrothermal scheduling problem  from
[Tutorial One: first steps](@ref) is:
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
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end
```

We will continue using this simple example of the hydrothermal scheduling problem. However, the formulation will be extended slightly.

We will add a terminal cost-to-go function to the formulation. The terminal cost-to-go is an important part of more developed hydrothermal scheduling problems. 

In the context of hydrothermal scheduling, the terminal cost is based on the concept of a marginal value of water. Water at the end-of-horizon has value because it can be used to generate electricity. The water is said to have a *marginal* value because a m^3 of water is worth more to when the reservoir is empty compared to when our reservoir is full.

The marginal water value function for this simple problem is shown in the chart below.

![marginal_water_value_chart](../assets/marginal_water_value_chart.png)

In the absence of a terminal cost-to-go function in the context of the hydrothermal scheduling problem, the policy would leave the reservoir empty at the end of the final stage. As water in hydro reservoirs has value, this outcome is undesirable and would not happen in reality. A terminal cost-to-go function is used to penalise such action. A terminal cost, which is a function of the final reservoir levels at the end of the final stage is included in the terminal stage objective. 

The terminal cost function, constructed from the marginal water value function used in the simple example is shown in the chart below:

![terminal_cost_chart](../assets/terminal_cost_chart.png)


## Formulating and solving the problem with SDDP

Now that the marginal water value and terminal cost function have been explained we can construct the model of the problem.

Compared to the simple hydrothermal scheduling problem presented in [Tutorial One: first steps](@ref) is the additional term in the objective in the final stage, the terminal cost-to-go.

```julia
using SDDP, JuMP, Clp

m = SDDPModel(
                 sense = :Min,
                stages = 3,
                solver = ClpSolver(),
       objective_bound = 0.0) do sp, t

    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        terminalcost       >= 0
     end)

    inflow = [50.0, 50.0, 50.0]
    fuel_cost = [50.0, 100.0, 150.0]
    terminal_marginal_cost = [-3 -2 -1]
    intercept = 600

    @constraints(sp, begin
        incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume
        thermal_generation + hydro_generation == 150
    end)
    for i in 1:length(terminal_marginal_cost)
        @constraint(sp, terminalcost >= terminal_marginal_cost[i] * outgoing_volume + intercept)
    end

    if t < 3
        @stageobjective(sp, fuel_cost[t] * thermal_generation)
    elseif t == 3
        @stageobjective(sp, fuel_cost[t] * thermal_generation + terminalcost)
    end
end

```
To solve this problem, we use the solve method:
```julia
status = solve(m; iteration_limit=5)
```

The output from the log is:
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
       15.600K         5.400K        |     1    0.0      0    0.0    0.0
        5.996K         5.600K        |     2    0.0      0    0.0    0.0
        5.600K         5.600K        |     3    0.0      0    0.0    0.0
        5.600K         5.600K        |     4    0.0      0    0.0    0.0
        5.600K         5.600K        |     5    0.0      0    0.0    0.0
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         5
        Termination Status: iteration_limit
===============================================================================
```

Including a terminal cost has increased the minimal policy cost from `5.0K` to `5.6K`. This additional cost of `0.6K` is due to the additional terminal cost term in the final stage objective.   


## Formulating the problem with infinite-horizon SDDP
In formulating many stochastic dynamic programs (such as the previous example), a terminating cost-to-go function is necessary. However, this terminating cost-to-go function is an assumption of many formulations, including the previous example. Solving a multi-stage stochastic dynamic problem with infinite-horizon stochastic dynamic programming (infinite-horizon SDDP) eliminates the need for a terminating cost-to-go function. 

The problem is constructed similarly to the problem in [Tutorial One: first steps](@ref). 
However, the first difference is in the input to [`SDDPModel`](@ref) method. The flag `is_infinite = true` tells [`SDDPModel`](@ref) to build the model using infinite-horizon SDDP. 

The additional inputs `lb_states` and `ub_states` provide the lower bound and upper bound on the state in the dummy stage 0 of the problem (explained next).

```julia
m = SDDPModel(
                 sense = :Min,
                stages = 3,
                solver = ClpSolver(),
       objective_bound = 0.0,
           is_infinite = true,
             lb_states = [0],
             ub_states = [200]) do sp, t
```

The next difference is containing all the standard constraints and the variables inside an `if` statement under the condition `if t > 0`. This is due to a dummy stage zero present under the hood. In the `else` statement the `@stage` should be constrained so the entering state to stage 1 (from the dummy stage zero), is as desired. In our simple example, in the dummy stage 0, we set the `outgoing_volume` to be equal to the `incoming_volume` leading to the incoming volume to stage 1 to be `200`. 

```julia
if t > 0
  # define all variables, , expressions, constraints, stage-objectives as in tutorial one
else
  @state(sp, outgoing_volume == incoming_volume, incoming_volume == 200)
end
```

Putting it all together, we get:

```julia
using SDDP, JuMP, Clp
m = SDDPModel(
                 sense = :Min,
                stages = 3,
                solver = ClpSolver(),
       objective_bound = 0.0,
           is_infinite = true,
             lb_states = [0],
             ub_states = [200]) do sp, t

    if t > 0
        @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
        @variables(sp, begin
            thermal_generation >= 0
            hydro_generation   >= 0
            hydro_spill        >= 0
            terminalcost       >= 0
         end)

        inflow = [50.0, 50.0, 50.0]
        fuel_cost = [50.0, 100.0, 150.0]

        @constraints(sp, begin
            incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume
            thermal_generation + hydro_generation == 150
        end)

        @stageobjective(sp, fuel_cost[t] * thermal_generation)
    else
        @state(sp, outgoing_volume == incoming_volume, incoming_volume == 200)
    end

end
```

To solve this problem, we use the [`solve`](@ref) method. The additional parameter `update_limit` is required to be passed to the solve function when using infinite-horizon SDDP. Choosing the values for `iteration_limit` and `update_limit` is more of an art than a science. 

For example, for a complex hydrothermal scheduling problem modelled with SDDP (with an exogenous terminal cost function), 3000 iterations of SDDP may be needed for convergence. When solving this problem with SDDP.jl, we would set `iteration_limit = 3000`. 

However, when the same problem is solved with infinite-horizon SDDP, a total of 8000 iterations of SDDP may be needed because the endogenous terminal cost-to-go function needs to converge. From my experience, we choose `iteration_limit = 500` and `update_limit = 20` (500x16 = 8000 total iteration of SDDP). Choosing the values for `iteration_limit` and `update limit` is discussed further in Section 5 of my [thesis](https://github.com/shasafoster/SDDP.jl/blob/master/docs/src/assets/foster_thesis.pdf). 

For a given problem usually solved with 

```julia
status = solve(m; iteration_limit=5, update_limit=10)
```

The output from the final update (update 10/10) log is:

```
-------------------------------------------------------------------------------
                          SDDP.jl © Oscar Dowson, 2017-2018
-------------------------------------------------------------------------------
    Solver:
        Serial solver
    Model:
        Stages:         4
        States:         1
        Subproblems:    4
        Value Function: Default
-------------------------------------------------------------------------------
              Objective              |  Cut  Passes    Simulations   Total
     Simulation       Bound   % Gap  |   #     Time     #    Time    Time
-------------------------------------------------------------------------------
        5.000K        17.500K        |     1    0.0      0    0.0    0.0
       22.500K        35.000K        |     2    0.0      0    0.0    0.0
       22.500K        35.000K        |     3    0.0      0    0.0    0.0
       22.500K        35.000K        |     4    0.0      0    0.0    0.0
       22.500K        35.000K        |     5    0.0      0    0.0    0.0
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         5
        Termination Status: iteration_limit
-------------------------------------------------------------------------------
Min Δ  : 12500.0
Max Δ  : 12500.0
Mean Δ : 12500.0
===============================================================================
```

Notice how the objective `Bound` is higher than the `Simulation` objective. This due to the when solving the problem with infinite-horizon SDDP overshoots the objective. The expected cost of the policy is more accurately shown through the value of Δ (the method of computing Δ is discussed in my [thesis](https://github.com/shasafoster/SDDP.jl/blob/master/docs/src/assets/foster_thesis.pdf)).

Note how the `Bound Objective - Δ = Simulation Objective`.

This concludes our tutorial 12 for SDDP.jl on infinite-horizon SDDP. 
 
## To do:

### Introduce uncertainty in dummy stage 0
There’s a minor issue with the algorithm. Currently there’s no uncertainty in the dummy stage 0. Because there is uncertatiny in all other stages, this causes issues when developing a proof demonstrating the infinite horizon SDDP methodology converges. The best solution would there to be no dummy stage 0 and have cuts generated from stage 1. 

### Hold cuts in memory
The current infinite-horizon methodology uses cuts written to a temporty directory (SDDP/src/temp/) vs holding the cuts in memory. Holding cuts in memory may result in faster solving, especially is `update_limit` is relatively large and `iteration_limit` is relatively small. 

### User-defined function to increase information value of the bound and simulation objective values
Given the `Simulation Objective` and the `Bound Objective` is less informative in the infinite-horizon programming, we output the Δ values as well. The values for Δ conveges as the policy converges. 

However, a user defined function that takes the initial `state` in the given iteration of SDDP as the input and applies the user-defined function to produce a string output would be useful in increasing the information value of the bound and simulation objective values procued from each iteration of SDDP.

For example for the hydrothermal sceduling problem, the `state` is commonly multi-dimensional. Let the state may be a vector of length 5 representing 5 reservoirs. The user defined function may multiply each of the of the vector entries by a specific power unique to each reservoir. This specific power converts the amount stored water in the each reservoir (m^3) to the amount stored hydroelectric energy (GWh) in each reservoir. The stored hydroelectric energies in each reservoir are summed to produce the net hydroelectric energy across the system. This value can then be written to the console with the `Simulation Objective` and the `Bound Objective` to inform these values. 

```julia
function write_informative_string(state::State)
    # do stuff with state
    return string
end
```

