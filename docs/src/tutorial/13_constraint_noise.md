# Tutorial Thirteen: constraint noise

In previous tutorials (e.g., [Tutorial Three: objective noise](@ref)), we
explicitly noted that SDDP.jl does not support noise in the constraint matrix.
However, we have recently developed a hack that works around this. It should not
be considered stable, and is an imperfect solution. The new version of JuMP will
enable native support for these modifications instead of the current hack.

!!!note
    This may break SDDP extensions such as [SDDiP.jl](https://github.com/lkapelevich/SDDiP.jl).
    It may also break some features of JuMP.

If `w ∈ [1,2,3]` with equal probability, and we want to add the constraint:
```julia
@constraint(sp, 2x + w*x <= 1)
```

Then in the subproblem definition, we can use the un-exported
[`SDDP.addconstraintnoise!`](@ref) method:
```julia
wx = SDDP.addconstraintnoise!(sp, x, [1,2,3])
@constraint(m, 2x + wx <= 1)
```
The first line introduces a new JuMP variable named `wx`. Depending on the
realization of the noise, this will take the value of `w*x`. Note that we named
the variable `wx`, but this decision was arbitrary.
[`SDDP.addconstraintnoise!`](@ref) is purposefully not exported to show that it
is an experimental feature.

!!!note
    This requires a solver that implements the [MathProgBase.changecoeffs!](http://mathprogbasejl.readthedocs.io/en/latest/lpqcqp.html#changecoeffs!)
    method.

To give a full example:
```julia
using SDDP, JuMP, Gurobi
m = SDDPModel(
        sense           = :Min,
        stages          = 2,
        objective_bound = 0.0,
        solver          = GurobiSolver(OutputFlag=0) ) do sp, t
    @state(sp, x>=0, x0==1.0)
    @stageobjective(sp, x)
    # instead of @constraint(sp, x == w * x0) where w ∈ [1,2], we write:
    r_x = SDDP.addconstraintnoise!(sp, x0, [1,2])
    @constraint(sp, x == r_x)
end
```
There are four possible scenarios that can be realized in this model:
 1. `(x₁, x₂) = (1, 1)` with total cost of `1+1=2`
 2. `(x₁, x₂) = (1, 2)` with total cost of `1+2=3`
 3. `(x₁, x₂) = (2, 2)` with total cost of `2+2=4`
 4. `(x₁, x₂) = (2, 4)` with total cost of `2+4=6`
 Therefore, the expected cost is `mean([2,3,4,6]) = 3.75`. If we solve the model
 for five iterations, we can confirm that the model converges to this bound:
```julia
julia> solve(m, iteration_limit=5)
 ... log omitted ...
julia> getbound(m)
 3.75
```
The log is:
```
-------------------------------------------------------------------------------
                          SDDP.jl © Oscar Dowson, 2017-2018
-------------------------------------------------------------------------------
    Solver:
        Serial solver
    Model:
        Stages:         2
        States:         1
        Subproblems:    2
        Value Function: Default
-------------------------------------------------------------------------------
              Objective              |  Cut  Passes    Simulations   Total
     Simulation       Bound   % Gap  |   #     Time     #    Time    Time
-------------------------------------------------------------------------------
        4.000          3.750         |     1    1.0      0    0.0    1.0
        3.000          3.750         |     2    1.0      0    0.0    1.3
        4.000          3.750         |     3    1.0      0    0.0    1.3
        2.000          3.750         |     4    1.0      0    0.0    1.3
        2.000          3.750         |     5    1.0      0    0.0    1.3
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         5
        Termination Status: iteration_limit
===============================================================================
```

Each uncertain term in the constraint matrix must be added using the
[`SDDP.addconstraintnoise!`](@ref) method. These terms can be used in
conjunction with right-hand side and  objective noise. In which case there must
be the same number of realizations in the constraint noise as there are in the
right-hand side noise terms. For example:
To give a full example:
```julia
using SDDP, JuMP, Gurobi
m = SDDPModel(
        sense           = :Min,
        stages          = 2,
        objective_bound = 0.0,
        solver          = GurobiSolver(OutputFlag=0) ) do sp, t
    @states(sp, begin
        x>=0, x0==1.0
        y>=0, y0==1.0
    end)
    # noise in the constraint matrix
    r_x = SDDP.addconstraintnoise!(sp, x0, [0.9,1.1])
    r_y = SDDP.addconstraintnoise!(sp, y0, [0.8,1.2])
    @constraint(sp, x == 0.9r_x + 0.1r_y)

    # noise in the right-hand side term
    @rhsnoise(sp, w=[0.0, 0.1], y == 0.9r_y + 0.1r_x + w)

    # noise in the objective function
    @stageobjective(sp, w=[1,2], w * x + y)
end
```

This concludes our thirteenth tutorial for SDDP.jl.
