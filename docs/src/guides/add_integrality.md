```@meta
CurrentModule = SDDP
```

# Integrality

There's nothing special about binary and integer variables in SDDP.jl. Your
models may contain a mix of binary, integer, or continuous state and control
variables. Use the standard JuMP syntax to add binary or integer variables.

For example:

```@example
using SDDP, HiGHS
model = SDDP.LinearPolicyGraph(
   stages = 3,
   lower_bound = 0.0,
   optimizer = HiGHS.Optimizer,
) do sp, t
   @variable(sp, 0 <= x <= 100, Int, SDDP.State, initial_value = 0)
   @variable(sp, 0 <= u <= 200, integer = true)
   @variable(sp, v >= 0)
   @constraint(sp, x.out == x.in + u + v - 150)
   @stageobjective(sp, 2u + 6v + x.out)
end
```

If you want finer control over how SDDP.jl computes subgradients in the backward
pass, you can pass an [`SDDP.AbstractDualityHandler`](@ref) to the
`duality_handler` argument of [`SDDP.train`](@ref).

See [Duality handlers](@ref) for the list of handlers you can pass.

## Convergence

SDDP.jl cannot guarantee that it will find a globally optimal policy when some
of the variables are discrete. However, in most cases we find that it can still
find an integer feasible policy that performs well in simulation.

Moreover, when the number of nodes in the graph is large, or there is
uncertainty, we are not aware of another algorithm that _can_ claim to find a
globally optimal policy.

## Does SDDP.jl implement the SDDiP algorithm?

Most discussions of SDDiP in the literature confuse two unrelated things.

 * First, how to compute dual variables
 * Second, when the algorithm will converge to a globally optimal policy
 
### Computing dual variables

The stochastic dual dynamic programming algorithm requires a subgradient of the
objective with respect to the incoming state variable. 

One way to obtain a valid subgradient is to compute an optimal value of the
dual variable ``\lambda`` in the following subproblem:

```math
\begin{aligned}
V_i(x, \omega) = \min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) + \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}[V_j(x^\prime, \varphi)]\\
& x^\prime = T_i(\bar{x}, u, \omega) \\
& u \in U_i(\bar{x}, \omega) \\
& \bar{x} = x \quad [\lambda],
\end{aligned}
```

The easiest option is to relax integrality of the discrete variables to form a
linear program and then use linear programming duality to obtain the dual. But
we could also use Lagrangian duality without needing to relax the integrality
constraints.

To compute the Lagrangian dual ``\lambda``, we penalize ``\lambda^\top(\bar{x} - x)``
in the objective instead of enforcing the constraint:
```math
\begin{aligned}
\max\limits_{\lambda}\min\limits_{\bar{x}, x^\prime, u} \;\; & C_i(\bar{x}, u, \omega) + \mathbb{E}_{j \in i^+, \varphi \in \Omega_j}[V_j(x^\prime, \varphi)] - \lambda^\top(\bar{x} - x)\\
& x^\prime = T_i(\bar{x}, u, \omega) \\
& u \in U_i(\bar{x}, \omega)
\end{aligned}
```

Compared with linear programming duality, the Lagrangian problem is difficult
to solve because it requires the solution of many mixed-integer programs
instead of a single linear program. This is one reason why "SDDiP" has poor
performance.

### Convergence

The second part to SDDiP is a very tightly scoped claim: _if_ all of the state
variables are binary _and_ the algorithm uses Lagrangian duality to compute a
subgradient, _then_ it will converge to an optimal policy.

In many cases, papers claim to "do SDDiP," but they have state variables which
are not binary. In these cases, the algorithm is not guaranteed to converge to a
globally optimal policy.

One work-around that has been suggested is to discretize the state variables
into a set of binary state variables. However, this leads to a large number of
binary state variables, which is another reason why "SDDiP" has poor
performance.

In general, we recommend that you introduce integer variables into your model
without fear of the consequences, and that you treat the resulting policy as a
good heuristic, rather than an attempt to find a globally optimal policy.
