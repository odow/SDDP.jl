# Tutorial Six: cut selection

Over the previous five tutorials, we formulated a simple risk-averse
hydrothermal scheduling problem. Now, in this tutorial, we improve the
computational performance of the solution process using *cut selection*.

As the SDDP algorithm progresses, many cuts (typically thousands) are added to
each subproblem. This increases the computational effort required to solve each
subproblem. In addition, many of the cuts created early in the solution process
may be redundant once additional cuts are added. This issue has spurred a vein
of research into heuristics for choosing cuts to keep or remove (this is
referred to as *cut selection*). To facilitate the development of cut selection
heuristics, SDDP.jl features the concept of a cut oracle associated with each
subproblem. A cut oracle has two jobs: it should store the complete list of cuts
created for that subproblem; and when asked, it should provide a subset of those
cuts to retain in the subproblem.

So far, only the Level One cut selection method is implemented. It can be
constructed using [`LevelOneCutOracle`](@ref). Cut selection can be added to a
model using the `cut_oracle` keyword argument to the [`SDDPModel`](@ref)
constructor.

Therefore, our risk-averse multistage stochastic optimization problem with cut
selection can be formulated in SDDP.jl as:
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

!!! info
    This will change once JuMP 0.19 lands.

Due to the design of JuMP, we are unable to delete cuts from the model.
Therefore, selecting a subset of cuts involves rebuilding the subproblems from
scratch. The user can control the frequency by which the cuts are selected and
the subproblems rebuilt with the `cut_selection_frequency` keyword argument of
the [`solve`](@ref) method. Frequent cut selection (i.e. when
`cut_selection_frequency` is small) reduces the size of the subproblems that are
solved but incurs the overhead of rebuilding the subproblems. However,
infrequent cut selection (i.e. when `cut_selection_frequency` is large) allows
the subproblems to grow large (by adding many constraints), leading to an
increase in the solve time of individual subproblems. Ultimately, this will be a
model-specific trade-off. As a rule of thumb, simpler (i.e. few variables and
constraints) models benefit from more frequent cut selection compared with
complicated (i.e. many variables and constraints) models.

```julia
status = solve(m,
    iteration_limit         = 10,
    cut_selection_frequency = 5
)
```

We can get the subproblem from a [`SDDPModel`](@ref) using
[`SDDP.getsubproblem`](@ref). For example, the subproblem in the first stage and
first Markov state is:
```julia
SDDP.getsubproblem(m, 1, 1)
```
Then, given a subproblem `sp`, we can get the cut oracle as follows:
```julia
oracle = SDDP.cutoracle(sp)
```
Finally, we can query the list of valid cuts in the oracle:
```julia
julia> SDDP.validcuts(oracle)
1-element Array{Cut, 1}:
 SDDP.Cut(29976.6, [-108.333])
```
So, despite performing 10 SDDP iterations, only one cut is needed to approximate
the cost-to-go function! A similar result can be found in the wet Markov state
in the second stage:
```julia
julia> SDDP.validcuts(SDDP.cutoracle(SDDP.getsubproblem(m, 2 1)))
3-element Array{Cut, 1}:
 SDDP.Cut(20625.0, [-162.5])
 SDDP.Cut(16875.0, [-112.5])
 SDDP.Cut(19375.0, [-137.5])
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

This concludes our sixth tutorial for SDDP.jl. In our next tutorial,
[Tutorial Seven: plotting](@ref), we discuss some of the plotting utilities of
SDDP.jl
