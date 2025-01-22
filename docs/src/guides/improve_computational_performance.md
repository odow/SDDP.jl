# Improve computational performance

SDDP is a computationally intensive algorithm. Here are some suggestions for
how the computational performance can be improved.

## Numerical stability (again)

We've already discussed this in the [Numerical stability](@ref) section of
[Words of warning](@ref). But, it's so important that we're going to
say it again: improving the problem scaling is one of the best ways to improve
the numerical performance of your models.

## Solver selection

The majority of the solution time is spent inside the low-level solvers.
Choosing which solver (and the associated settings) correctly can lead to big
speed-ups.

 - Choose a commercial solver.

   Options include [CPLEX](https://github.com/jump-dev/CPLEX.jl),
   [Gurobi](https://github.com/jump-dev/Gurobi.jl), and
   [Xpress](https://github.com/jump-dev/Xpress.jl). Using free solvers such as
   [CLP](https://github.com/jump-dev/Clp.jl) and
   [HiGHS](https://github.com/jump-dev/HiGHS.jl) isn't a viable approach for large
   problems.

- Try different solvers.

  Even commercial solvers can have wildly different solution times. We've seen
  models on which CPLEX was 50% fast than Gurobi, and vice versa.

- Experiment with different solver options.

  Using the default settings is usually a good option. However, sometimes it can
  pay to change these. In particular, forcing solvers to use the dual simplex
  algorithm (e.g., [`Method=1` in Gurobi](https://www.gurobi.com/documentation/8.1/refman/method.html)
  ) is usually a performance win.

## Single-cut vs. multi-cut

There are two competing ways that cuts can be created in SDDP: _single_-cut and
_multi_-cut. By default, `SDDP.jl` uses the _single-cut_ version of SDDP.

The performance of each method is problem-dependent. We recommend that you try
both in order to see which one performs better. In general, the _single_-cut
method works better when the number of realizations of the stagewise-independent
random variable is large, whereas the multi-cut method works better on small
problems. However, the multi-cut method can cause numerical stability problems,
particularly if used in conjunction with objective or belief state variables.

You can switch between the methods by passing the relevant flag to `cut_type` in
[`SDDP.train`](@ref).
```julia
SDDP.train(model; cut_type = SDDP.SINGLE_CUT)
SDDP.train(model; cut_type = SDDP.MULTI_CUT)
```

## Parallelism

SDDP.jl can take advantage of the parallel nature of modern computers to solve
problems across multiple threads.

Start Julia from a command line with `N` threads using the `--threads` flag:
```julia
julia --threads N
```

Then, pass an instance of [`SDDP.Threaded`](@ref) to the `parallel_scheme`
argument of [`SDDP.train`](@ref) and [`SDDP.simulate`](@ref).

```julia
using SDDP, HiGHS
model = SDDP.LinearPolicyGraph(
  stages = 2, lower_bound = 0, optimizer = HiGHS.Optimizer
) do sp, t
     @variable(sp, x >= 0, SDDP.State, initial_value = 1)
     @stageobjective(sp, x.out)
end
SDDP.train(model; iteration_limit = 10, parallel_scheme = SDDP.Threaded())
SDDP.simulate(model, 10; parallel_scheme = SDDP.Threaded())
```
