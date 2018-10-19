# Tutorial 15: Markovian policy graphs with infinite-horizon SDDP

Previously in tutorial [Tutorial Fourteen: Infinite Horizon](@ref) we discussed how to formulate a problem with infinite-horizon stochastic dual dynamic programming. When using infinite-horizon SDDP introducing *stagewise-dependent* uncertainty using a Markov chain follows a slightly different method than in [Tutorial Three: objective noise](@ref). This is because when using infinite-horizon SDDP we are finding the steady state policy.

Recall that our model for the hydrothermal scheduling problem  from
[Tutorial Fourteen: Infinite Horizon](@ref) is:
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
             
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
    )
    setnoiseprobability!(sp, [1/3, 1/3, 1/3])
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end
```

## Formulating the problem

As in [Tutorial Four: Markovian policy graphs](@ref)), we consider a Markov chain 
with two *climate* states: wet and dry.

The difference between the Markov method with standard SDDP and Markov methods with infinite-horizon SDDP is that there is a  transition matrix between stage T and stage 1 in infinite-horizon SDDP instead on a transition matrix between stage '0' and 1 and in standard SDDP as in [Tutorial Four: Markovian policy graphs](@ref). 

This is because in [Tutorial Four: Markovian policy graphs](@ref) we start from the same Markov state (state 1) in each iteration of SDDP. However in infinite-horizon SDDP, we the conclusion of one iteration leads into the start of the next as the next iteration will start from where the previous iteration finished. Precisely, the inital state and markov state in stage 1 of  iteration *i* is the same as the final state and markov state at the end of stage T in iteration *i-1*. 

Hence, the vector of
Markov transition matrices is given by:
```julia
T = Array{Float64, 2}[
    [ 0.75 0.25 ; 0.25 0.75 ],
    [ 0.75 0.25 ; 0.25 0.75 ],
    [ 0.75 0.25 ; 0.25 0.75 ]
]
```

To add the Markov chain to the model, modifications are required. First, we
give the vector of transition matrices to the [`SDDPModel`](@ref) constructor
using the `markov_transition` keyword. Second, the `do sp, t ... end` syntax
is extended to `do sp, t, i ... end`, where `i` is the index of the Markov state
and runs from `i=1` to the number of Markov states in stage `t`. Now, both `t`
and `i` can be used anywhere inside the subproblem definition.

Our model is now:
```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0,
      markov_transition = Array{Float64, 2}[
          [ 1.0 ]',
          [ 0.75 0.25 ],
          [ 0.75 0.25 ; 0.25 0.75 ]
      ]
                                        ) do sp, t, i
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
    )
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],
        mupliplier * fuel_cost[t] * thermal_generation
    )
    if i == 1  # wet climate state
        setnoiseprobability!(sp, [1/6, 1/3, 0.5])
    else       # dry climate state
        setnoiseprobability!(sp, [0.5, 1/3, 1/6])
    end
end
```

## Solving the problem

Now we need to solve the problem. As in the previous two tutorials, we use the
[`solve`](@ref) function. However, this time we terminate the SDDP algorithm
by setting a time limit (in seconds) using the `time_limit` keyword:
```julia
status = solve(m,
    time_limit = 0.05
)
```
The termination status is `:time_limit`, and the output from the log is
now:
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
        0.000          6.198K        |     1    0.0      0    0.0    0.0
        2.000K         7.050K        |     2    0.0      0    0.0    0.0
        2.000K         7.050K        |     3    0.0      0    0.0    0.0
        2.000K         7.135K        |     4    0.0      0    0.0    0.0
        5.000K         7.135K        |     5    0.0      0    0.0    0.0
        2.000K         7.135K        |     6    0.0      0    0.0    0.0
        2.000K         7.135K        |     7    0.0      0    0.0    0.0
        2.000K         7.135K        |     8    0.0      0    0.0    0.0
        2.000K         7.135K        |     9    0.0      0    0.0    0.0
        9.000K         7.135K        |    10    0.0      0    0.0    0.0
        2.000K         7.135K        |    11    0.0      0    0.0    0.0
        5.000K         7.135K        |    12    0.1      0    0.0    0.1
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         12
        Termination Status: time_limit
===============================================================================
```

## Understanding the solution

Instead of performing a Monte Carlo simulation, you may want to simulate one
particular sequence of noise realizations. This *historical* simulation can
also be conducted using the [`simulate`](@ref) function.
```julia
simulation_result = simulate(m,
    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill],
    noises       = [1, 1, 3],
    markovstates = [1, 2, 2]
)
```
Again, `simulation_result` is a single dictionary. In addition to the variable
values and the special keys `:noise`, `:objective`, and `:stageobjective`,
SDDP.jl also records the index of the Markov state in each stage via the
`:markov` key. We can confirm that the historical sequence of Markov states was
visited as follows:
```julia
julia> simulation_result[:markov]
3-element Array{Int, 1}:
 1
 2
 2
```

This concludes our fourth tutorial for SDDP.jl. In the next tutorial,
[Tutorial Five: risk](@ref), we introduce risk into the problem.
