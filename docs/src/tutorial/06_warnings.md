# Basics VI: words of warning

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
```jldoctest; filter=r"\|.+?\n"
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

Kokako.train(model, iteration_limit = 5)

println("Finished training!")

# output

———————————————————————————————————————————————————————————————————————————————
                        Kokako - © Oscar Dowson, 2018-19.
———————————————————————————————————————————————————————————————————————————————
 Iteration | Simulation |      Bound |   Time (s)
———————————————————————————————————————————————————————————————————————————————
         1 |     6.500  |     3.000  |     0.001
         2 |     3.500  |     3.500  |     0.001
         3 |     3.500  |     3.500  |     0.004
         4 |     3.500  |     3.500  |     0.004
         5 |     3.500  |     3.500  |     0.005
Finished training!
```

Now consider the case when we set the `lower_bound` to `10.0`:

```jldoctest; filter=r"\|.+?\n"
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

Kokako.train(model, iteration_limit = 5)

println("Finished training!")

# output

———————————————————————————————————————————————————————————————————————————————
                        Kokako - © Oscar Dowson, 2018-19.
———————————————————————————————————————————————————————————————————————————————
 Iteration | Simulation |      Bound |   Time (s)
———————————————————————————————————————————————————————————————————————————————
         1 |     6.500  |    11.000  |     0.001
         2 |     5.500  |    11.000  |     0.002
         3 |     5.500  |    11.000  |     0.002
         4 |     5.500  |    11.000  |     0.003
         5 |     5.500  |    11.000  |     0.006
Finished training!
```

How do we tell which is more appropriate? There are a few clues that you should
look out for.

- The bound converges to a value above (if minimizing) the simulated cost of the
  policy. In this case, the problem is deterministic, so it is easy to tell. But
  you can also check by performing a Monte Carlo simulation like we did in
  [Basics II: adding uncertainty](@ref).

- The bound converges to different values when we change the bound. This is
  another clear give-away. The bound provided by the user is only used in the
  initial iterations. __It should not change the value of the converged
  policy.__ Thus, if you don't know an appropriate value for the bound, choose
  an initial value, and then increase (or decrease) the value of the bound to
  confirm that the value of the policy doesn't change.

- The bound converges to a value _close_ to the bound provided by the user. This
  varies between models, but notice that `11.0` is quite close to `10.0`
  compared with `3.5` and `0.0`.

This concludes or series of basic introductory tutorials for `Kokako.jl`. When
you're ready, continue to our intermediate series of tutorials, beginning with
[Intermediate I: risk](@ref).
