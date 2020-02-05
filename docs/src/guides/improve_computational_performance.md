# Improve computational performance

SDDP is a computationally intensive algorithm. Here are some suggestions for
how the computational performance can be improved.

## Numerical stability (again)

We've already discussed this in the [Numerical stability](@ref) section of
[Basic VI: words of warning](@ref). But, it's so important that we're going to
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

### Single-cut vs. multi-cut

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

### Parallelism

SDDP.jl can take advantage of the parallel nature of modern computers to solve problems
across multiple cores.

!!! info
    We highly recommend that you read the Julia manual's section on [parallel computing](https://docs.julialang.org/en/v1/manual/parallel-computing/).

You can start Julia from a command line with `N` processors using the `-p` flag:
```julia
julia -p N
```

Alternatively, you can use the `Distributed.jl` package:
```julia
using Distributed
Distributed.addprocs(N)
```

!!! warning
    Workers **DON'T** inherit their parent's Pkg environment. Therefore, if you started
    Julia with `--project=/path/to/environment` (or if you activated an environment from the
    REPL), you will need to put the following at the top of your script:
    ```julia
    using Distributed
    @everywhere begin
        import Pkg
        Pkg.activate("/path/to/environment")
    end
    ```

Currently SDDP.jl supports to parallel schemes, [`SDDP.Serial`](@ref) and
[`SDDP.Asynchronous`](@ref). Instances of these parallel schemes should be passed to the
`parallel_scheme` argument of [`SDDP.train`](@ref) and [`SDDP.simulate`](@ref).

```julia
using SDDP, GLPK
model = SDDP.LinearPolicyGraph(
  stages = 2, lower_bound = 0, optimizer = GLPK.Optimizer
) do sp, t
     @variable(sp, x >= 0, SDDP.State, initial_value = 1)
     @stageobjective(sp, x.out)
end
SDDP.train(model; iteration_limit = 10, parallel_scheme = SDDP.Asynchronous())
SDDP.simulate(model, 10; parallel_scheme = SDDP.Asynchronous())
```

There is a large overhead for using the asynchronous solver. Even if you choose asynchronous
mode, SDDP.jl will start in serial mode while the initialization takes place. Therefore, in
the log you will see that the initial iterations take place on the master thread (`Proc. ID
= 1`), and it is only after while that the solve switches to full parallelism.

!!! info
    Because of the large data communication requirements (all cuts have to be shared with
    all other cores), the solution time will not scale linearly with the number of cores.

!!! info
    Given the same number of iterations, the policy obtained from asynchronous mode will be
    _worse_ than the policy obtained from serial mode. However, the asynchronous solver can
    take significantly less time to compute the same number of iterations.

#### Data movement

By defualt, data defined on the master process is not made available to the workers.
Therefore, a model like the following:
```julia
data = 1
model = SDDP.LinearPolicyGraph(stages = 2, lower_bound = 0) do sp, t
     @variable(sp, x >= 0, SDDP.State, initial_value = data)
     @stageobjective(sp, x.out)
end
```
will result in an `UndefVarError` error like `UndefVarError: data not defined`.

There are three solutions for this problem.

##### Option 1: declare data inside the build function

```julia
model = SDDP.LinearPolicyGraph(stages = 2) do sp, t
    data = 1
    @variable(sp, x >= 0, SDDP.State, initial_value = 1)
    @stageobjective(sp, x)
end
```

##### Option 2: use `@everywhere`

```julia
@everywhere begin
    data = 1
end
model = SDDP.LinearPolicyGraph(stages = 2) do sp, t
    @variable(sp, x >= 0, SDDP.State, initial_value = 1)
    @stageobjective(sp, x)
end
```

##### Option 3: build the model in a function

```julia
function build_model()
    data = 1
    return SDDP.LinearPolicyGraph(stages = 2) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 1)
        @stageobjective(sp, x)
    end
end

model = build_model()
```

#### Initialization hooks

[`SDDP.Asynchronous`](@ref) accepts a pre-processing hook that is run on each worker process
_before_ the model is solved. The most useful situation is for solvers than need an
initialization step. A good example is Gurobi, which can share an environment amongst all
models on a worker. Notably, this environment **cannot** be shared amongst workers, so
defining one environment at the top of a script will fail!

To initialize a new environment on each worker, use the following:
```julia
SDDP.train(
    model,
    parallel_scheme = SDDP.Asynchronous() do m
        env = Gurobi.Env()
        optimizer = with_optimizer(Gurobi.Optimizer, env, OutputFlag = 0)
        for node in values(m.nodes)
            set_optimizer(node.subproblem, optimizer)
        end
    end
)
```
