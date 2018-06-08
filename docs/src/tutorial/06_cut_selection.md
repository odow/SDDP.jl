# Tutorial Six: cut selection

Over the previous five tutorials, we formulated a simple risk-averse
hydrothermal scheduling problem. Now, in this tutorial, we improve the
computational performance of the solution process using *cut selection*.

Recall that our model for the hydrothermal scheduling problem from
[Tutorial Five: risk](@ref) is:

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
      ],
           risk_measure = EAVaR(lambda=0.5, beta=0.1)
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

## Formulating the problem




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
      ],
           risk_measure = EAVaR(lambda=0.5, beta=0.1),
             cut_oracle = LevelOneCutOracle()
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

```julia
status = solve(m,
    max_iterations          = 10,
    cut_selection_frequency = 5
)
```

```
-------------------------------------------------------------------------------
                          SDDP.jl © Oscar Dowson, 2017-2018
-------------------------------------------------------------------------------
    Solver:
        Serial solver
    Model:
        Stages:         3
        States:         1
        Subproblems:    5
        Value Function: Default
-------------------------------------------------------------------------------
              Objective              |  Cut  Passes    Simulations   Total
     Simulation       Bound   % Gap  |   #     Time     #    Time    Time
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         10
        Termination Status: max_iterations
===============================================================================
```

## Extra for experts: new cut selection heurisitics

One of the cool features of SDDP.jl is how easy it is to create new cut
selection heuristics.

To illustrate this, we implement the *Last Cuts* strategy. This heuristic
selects the last `N` discovered cuts to keep in each subproblem. First, we need
to create a new concrete subtype of the abstract type `AbstractCutOracle`
defined by SDDP.jl:
```julia
"""
    LastCuts(N::Int)

Create a cut oracle that keeps the last `N` discovered cuts.
"""
struct LastCuts <: SDDP.AbstractCutOracle
    cuts::Vector{SDDP.Cut}
    N::Int
    LastCuts(N::Int) = new(SDDP.Cut[], N)
end
```

Then, we need to overload two methods: `SDDP.storecut!` and `SDDP.validcuts`.

`SDDP.storecut` takes four arguments. The first is an instance of the cut
oracle. The second and third arguments are the SDDPModel `m` and the JuMP
subproblem `sp`. The fourth argument is the cut itself. For our example, we
store all the cuts that have been discovered:
```julia
function SDDP.storecut!(oracle::LastCuts, m::SDDPModel, sp::JuMP.Model, cut::SDDP.Cut)
    push!(oracle.cuts, cut)
end
```

`SDDP.validcuts` returns a list of all of the cuts that are valid cuts to keep
in the subproblem. The `LastCuts` oracle returns the most recent `N` discovered
cuts:
```julia
function SDDP.validcuts(oracle::LastCuts)
    oracle.cuts[max(1, end - oracle.N + 1):end]
end
```

Now, the cut oracle `LastCuts(500)` can be used just like any other cut oracle
defined by SDDP.jl.

This concludes our sixth tutorial for SDDP.jl.
