# Intermediate V: performance

SDDP is a computationally intensive algorithm. In this tutorial, we give
suggestions for how the computational performance can be improved.

## Numerical stability (again)

We've already discussed this in the [Numerical stability](@ref) section of
[Basics VI: words of warning](@ref). But, it's so important that we're going to
say it again: improving the problem scaling is one of the best ways to improve
the numerical performance of your models.

## Solver selection

The majority of the solution time is spent inside the low-level solvers.
Choosing which solver (and the associated settings) correctly can lead to big
speed-ups.

 - Choose a commercial solver.

   Options include [CPLEX](https://github.com/JuliaOpt/CPLEX.jl),
   [Gurobi](https://github.com/JuliaOpt/Gurobi.jl), and
   [Xpress](https://github.com/JuliaOpt/Xpress.jl). Using free solvers such as
   [CLP](https://github.com/JuliaOpt/Clp.jl) and
   [GLPK](https://github.com/JuliaOpt/GLPK.jl) isn't a viable approach for large
   problems.

- Try different solvers.

  Even commercial solvers can have wildly different solution times. We've seen
  models on which CPLEX was 50% fast than Gurobi, and vice versa.

- Experiment with different solver options.

  Using the default settings is usually a good option. However, sometimes it can
  pay to change these. In particular, forcing solvers to use the dual simplex
  algorithm (e.g., [`Method=1` in Gurobi](https://www.gurobi.com/documentation/8.1/refman/method.html)
  ) is usually a performance win.
