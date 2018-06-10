# Tutorial Nine: nonlinear models

In our previous tutorials, we formulated a linear version of the hydrothermal
scheduling problem. To do so, we had to make a large assumption, namely, that
the inflows were *stagewise-independent*. In
[Tutorial Four: Markovian policy graphs](@ref), we improved upon this slightly
using a Markov chain with persistent *climate* states. However, another
commonly used model is to assume that the inflows follow a log auto-regressive
process. In this tutorial, we explain how to implement this in SDDP.jl.

As a starting point, we use a model that is very similar to the hydrothermal
scheduling problem we formulated in [Tutorial Two: RHS noise](@ref). In that
model, we assumed that there the a constraint:
```julia
@rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
    outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
)
```
We assumed that the `inflow` term was *stagewise-independent*. In this tutorial,
we assume that the inflow term can be modelled like:
```
log(inflow[t]) = log(inflow[t-1]) + log(noise),
```
where `noise` is drawn from `[0.9, 1.0, 1.1]` with uniform probability. Since we
value of the inflow in the previous stage affects the value in the current
stage, it needs to be added as a state variable:
```julia
@state(sp, current_inflow >= 1, previous_inflow == 50)
```
Note that we place a small, positive lower bound on the `current_inflow` to
prevent the solver calling the undefined `log(0.0)`.

Now we want to add the log transition constraint. We can do this using JuMP's
`@NLconstraint` macro to add this constraint. However, SDDP.jl only supports
noise terms in the right hand-side of a linear constraint. We can overcome this
limitation by introducing a dummy variable (note that it also has a small,
positive lower bound to avoid `log(0.0)`):
```julia
@variable(sp, noise >= 0.1)
@rhsnoise(sp, ω = [0.9, 1.0, 1.1], noise == ω)
@NLconstraint(sp, log(current_inflow) == log(previous_inflow) + log(noise))
```

Finally, we're using JuMP's nonlinear functionality, we need to choose an
appropriate solver. We choose to use the COIN solver
[Ipopt](https://github.com/JuliaOpt/Ippopt.jl). Therefore, our final model is:
```julia
using JuMP, SDDP, Ipopt
m = SDDPModel(
    stages          = 3,
    sense           = :Min,
    solver          = IpoptSolver(print_level=0),
    objective_bound = 0.0
                ) do sp, t
    @states(sp, begin
        200 >= outgoing_volume >= 0, incoming_volume == 200
               current_inflow  >= 1, previous_inflow ==  50
    end)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        inflow_noise_term  >= 0.1
    end)
    @constraints(sp, begin
        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == current_inflow
        hydro_generation + thermal_generation >= 150.0
    end)

    @rhsnoise(sp, ω = [0.9, 1.0, 1.1], inflow_noise_term == ω)
    @NLconstraint(sp, log(current_inflow) == log(previous_inflow) + log(inflow_noise_term))

    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end
```

This problem can be solved just like any other SDDP model:
```julia
status = solve(m, iteration_limit = 10)
simulation_result = simulate(m, 1, [:inflow′])
```
Then, we can check that the inflows do indeed follow a log auto-regressive
process.
```julia
julia> simulaton_result[1][:inflow′]
 55.0
 60.5
 54.45
```

That concludes our ninth tutorial for SDDP.jl. In our next tutorial,
[Tutorial Ten: parallelism](@ref), we explain how to solve SDDP models in parallel.
