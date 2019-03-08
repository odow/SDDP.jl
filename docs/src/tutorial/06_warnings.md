```@meta
CurrentModule = Kokako
```

# Basic VI: words of warning

SDDP is a powerful solution technique for multistage stochastic programming.
However, there are a number of subtle things to be aware of before creating
your own models.

## Numerical stability

If you aren't aware, SDDP builds an outer-approximation to a convex function
using cutting planes. This results in a formulation that is particularly hard
for solvers like GLPK, Gurobi, and CPLEX to deal with. As a result, you may run
into weird behavior. This behavior could include

 - Iterations suddenly taking a long time (the solver stalled)
 - Subproblems turning infeasible or unbounded after many iterations
 - Solvers returning "Numerical Error" statuses

### Problem scaling

In almost all cases, the cause of this is poor problem scaling. For our purpose,
poor problem scaling means having variables with very large numbers and
variables with very small numbers in the same model.

!!! info
    Gurobi has an excellent [set of articles](http://www.gurobi.com/documentation/8.1/refman/numerics_gurobi_guidelines.html)
    on numerical issues and how to avoid them.

Consider, for example, the hydro-thermal scheduling problem we have been
discussing in previous tutorials.

If we define the volume of the reservoir in terms of m³, then a lake might have
a capacity of 10^10 m³: `@variable(subproblem, 0 <= volume <= 10^10)`. Moreover,
the cost per cubic meter might be around \\\$0.05/m³. To calculate the  value of
water in our reservoir, we need to multiple a variable on the order of 10^10, by
one on the order of 10⁻²!. That is twelve orders of magnitude!

To improve the performance of the SDDP algorithm (and reduce the chance of weird
behavior), try to re-scale the units of the problem in order to reduce the
largest difference in magnitude. For example, if we talk in terms of million m³,
then we have a capacity of 10⁴ million m³, and a price of \\\$50,000 per million
m³. Now things are only one order of magnitude apart.

### Numerical stability report

To aid in the diagnose of numerical issues, you can call
[`Kokako.numerical_stability_report`](@ref). By default, this aggregates all of
the nodes into a single report. You can produce a stability report for each node
by passing `by_node=true`.

```jldoctest
using Kokako

model = Kokako.LinearPolicyGraph(
        stages = 2, lower_bound = -1e10, direct_mode=false) do subproblem, t
    @variable(subproblem, x >= -1e7, Kokako.State, initial_value=1e-5)
    @constraint(subproblem, 1e9 * x.out >= 1e-6 * x.in + 1e-8)
    @stageobjective(subproblem, 1e9 * x.out)
end

Kokako.numerical_stability_report(model)

# output

Numerical stability report
  Non-zero Matrix range     [1e-06, 1e+09]
  Non-zero Objective range  [1e+00, 1e+09]
  Non-zero Bounds range     [1e+07, 1e+10]
  Non-zero RHS range        [1e-08, 1e-08]
WARNING: numerical stability issues detected
  - Matrix range contains small coefficients
  - Matrix range contains large coefficients
  - Objective range contains large coefficients
  - Bounds range contains large coefficients
  - RHS range contains small coefficients
These coefficients can cause numerical stability issues. Consider reformulating
the model.
```

By default, a numerical stability check is run when you call
[`Kokako.train`](@ref), although it can be turned off by passing
`perform_numerical_stability_check = false`.

### Solver-specific options

If you have a particularly troublesome model, you should investigate setting
solver-specific options to improve the numerical stability of each solver. For
example, Gurobi has a [`NumericFocus` option](http://www.gurobi.com/documentation/8.1/refman/numericfocus.html#parameter:NumericFocus).

## Choosing an initial bound

One of the important requirements when building a SDDP model is to choose an
appropriate bound on the objective (lower if minimizing, upper if maximizing).
However, it can be hard to choose a bound if you don't know the solution! (Which
is very likely.)

The bound should not be as large as possible (since this will help with
convergence and the numerical issues discussed above), but if chosen to small,
it may cut of the feasible region and lead to a sub-optimal solution.

Consider the following simple model, where we first set `lower_bound` to `0.0`.
```jldoctest
using Kokako, GLPK

model = Kokako.LinearPolicyGraph(
            stages = 3,
            sense = :Min,
            lower_bound = 0.0,
            optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, t
    @variable(subproblem, x >= 0, Kokako.State, initial_value = 2)
    @variable(subproblem, u >= 0)
    @variable(subproblem, v >= 0)
    @constraint(subproblem, x.out == x.in - u)
    @constraint(subproblem, u + v == 1.5)
    @stageobjective(subproblem, t * v)
end

Kokako.train(model, iteration_limit = 5, perform_numerical_stability_check=false)

# output

----------------------------------------------------
         SDDP.jl (c) Oscar Dowson, 2017-19

 Iteration   Simulation      Bound        Time (s)
         1   6.50000e+00   3.00000e+00   0.00000e+00
         2   3.50000e+00   3.50000e+00   1.00017e-03
         3   3.50000e+00   3.50000e+00   1.00017e-03
         4   3.50000e+00   3.50000e+00   2.00009e-03
         5   3.50000e+00   3.50000e+00   2.00009e-03

Terminating training with status: iteration_limit
----------------------------------------------------
```

Now consider the case when we set the `lower_bound` to `10.0`:

```jldoctest
using Kokako, GLPK

model = Kokako.LinearPolicyGraph(
            stages = 3,
            sense = :Min,
            lower_bound = 10.0,
            optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, t
    @variable(subproblem, x >= 0, Kokako.State, initial_value = 2)
    @variable(subproblem, u >= 0)
    @variable(subproblem, v >= 0)
    @constraint(subproblem, x.out == x.in - u)
    @constraint(subproblem, u + v == 1.5)
    @stageobjective(subproblem, t * v)
end

Kokako.train(model, iteration_limit = 5, perform_numerical_stability_check=false)

# output

----------------------------------------------------
         SDDP.jl (c) Oscar Dowson, 2017-19

 Iteration   Simulation      Bound        Time (s)
         1   6.50000e+00   1.10000e+01   0.00000e+00
         2   5.50000e+00   1.10000e+01   0.00000e+00
         3   5.50000e+00   1.10000e+01   9.99928e-04
         4   5.50000e+00   1.10000e+01   9.99928e-04
         5   5.50000e+00   1.10000e+01   9.99928e-04

Terminating training with status: iteration_limit
----------------------------------------------------
```

How do we tell which is more appropriate? There are a few clues that you should
look out for.

- The bound converges to a value above (if minimizing) the simulated cost of the
  policy. In this case, the problem is deterministic, so it is easy to tell. But
  you can also check by performing a Monte Carlo simulation like we did in
  [Basic II: adding uncertainty](@ref).

- The bound converges to different values when we change the bound. This is
  another clear give-away. The bound provided by the user is only used in the
  initial iterations. __It should not change the value of the converged
  policy.__ Thus, if you don't know an appropriate value for the bound, choose
  an initial value, and then increase (or decrease) the value of the bound to
  confirm that the value of the policy doesn't change.

- The bound converges to a value _close_ to the bound provided by the user. This
  varies between models, but notice that `11.0` is quite close to `10.0`
  compared with `3.5` and `0.0`.

This concludes our sixth tutorial for `SDDP.jl`. In our final basics tutorial,
[Basic VII: modelling tips](@ref) we discuss some additional modelling tips.
