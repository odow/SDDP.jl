# Markovian policy graphs

In our three tutorials ([First steps](@ref), [RHS noise](@ref), and [Objective
noise](@ref)), we formulated a simple hydrothermal scheduling problem with
stagewise-independent noise in the right-hand side of the constraints and in the
objective function. Now, in this tutorial, we introduce some
_stagewise-dependent_ uncertainty using a Markov chain.

Recall that our model for the hydrothermal scheduling problem  from
[RHS noise](@ref) is:
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
    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
    )
    setnoiseprobability!(sp, [1/3, 1/3, 1/3])
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],
        mupliplier * fuel_cost[t] * thermal_generation
    )
end
```

## Formulating the problem

In this tutorial we consider a Markov chain with two _climate_ states: wet and
dry. Each Markov state is associated with an integer, in this case the wet
climate state  is Markov state `1` and the dry climate state is Markov state
`2`. In the wet climate state, the probability of the high inflow increases to
50%, and the probability of the low inflow decreases to 1/6. In the dry climate
state, the converse happens. There is also persistence in the climate state: the
probability of remaining in the current state is 75%, and the probability of
transitioning to the other climate state is 25%. We assume that the first stage
starts in the _wet_ climate state.

For each stage, we need to provide a Markov transition matrix. This is an
`M`x`N` matrix, where the element `A[i,j]` gives the probability of transitioning
from Markov state `i` in the previous stage to Markov state `j` in the current
stage. The first stage is special because we assume there is a "zero'th" stage
which has one Markov state. Furthermore, the number of columns in the transition
matrix of a stage (i.e., the number of Markov states) must equal the number of
rows in the next stage's transition matrix. For our example, the vector of
Markov transition matrices is given by:
```julia
T = Array{Float64, 2}[
    [ 1.0 0.0 ],
    [ 0.75 0.25 ; 0.25 0.75 ],
    [ 0.75 0.25 ; 0.25 0.75 ]
]
```
However, note that we never sample the dry Markov state in stage one. Therefore,
we can drop that Markov state so that there is only one Markov state in stage 1.
We also need to modify the transition matrix in stage 2 to account for this:
```julia
T = Array{Float64, 2}[
    [ 1.0 ]',
    [ 0.75 0.25 ],
    [ 0.75 0.25 ; 0.25 0.75 ]
]
```

To add the Markov chain to the model, we modifications are required. First, we
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
particular sequence of noise realizations. This _historical_ simulation can
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
[Risk](@ref), we introduce risk into the problem.
