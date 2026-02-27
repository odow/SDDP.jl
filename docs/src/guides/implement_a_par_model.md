# How to implement a PAR model

A common stochastic model is the Periodic Autoregressive model:
```math
y_t = \sum\limits_{p=1}^P \alpha_p y_{t-p} + \omega_t
```
To implement this model, you must first:

 * compute the constants $\alpha_p$ separately
 * construct a distribution for the noise term $\omega_t$

Then, implement the SDDP model as follows:

```@repl
using SDDP, HiGHS
alpha = [0.5, 0.3, 0.2]
y_history = [1.0, 1.0, 1.0]
Ω = [-0.25, 0.0, 0.25]
model = SDDP.LinearPolicyGraph(;
    stages = 12,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    P = length(alpha)
    @variable(sp, y[1:P], SDDP.State, initial_value = 1.0)
    @variable(sp, omega)
    # Add the constraint for the PAR model
    @constraint(sp, y[1].out == sum(alpha[p] * y[p].in for p in 1:P) + omega)
    # Add a "shift" constraint to move the values of `y` between stages
    @constraint(sp, [p in 2:P], y[p].out == y[p-1].in)
    # Implement the PAR model by parameterizing the additive noise term
    SDDP.parameterize(sp, Ω) do ω
        fix(omega, ω; force = true)
    end
end
simulations = SDDP.simulate(model, 1, [:y]);
y = [stage[:y][1].out for stage in simulations[1]]
sampling_scheme = SDDP.Historical(tuple.(1:3, [0.1, -0.1, 0.2]));
simulations = SDDP.simulate(model, 1, [:y]; sampling_scheme);
y = [stage[:y][1].out for stage in simulations[1]]
```

## Use a more sophisticated noise term that allows simulating actual sequences

If you want the ability to sample additive noise terms during training, but
actual historical sequences during simulation, then you can use a more
sophisticated noise term that allows you to control the `SDDP.parameterize`
function:

```@repl
using SDDP, HiGHS
alpha = [0.5, 0.3, 0.2]
y_history = [1.0, 1.0, 1.0]
struct NoiseTerm
    is_noise::Bool
    value::Float64
end
Ω = NoiseTerm.(true, [-0.25, 0.0, 0.25])
model = SDDP.LinearPolicyGraph(;
    stages = 12,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    P = length(alpha)
    @variable(sp, y[1:P], SDDP.State, initial_value = 1.0)
    @variable(sp, omega)
    @constraint(sp, y[1].out == sum(alpha[p] * y[p].in for p in 1:P) + omega)
    @constraint(sp, [p in 2:P], y[p].out == y[p-1].in)
    SDDP.parameterize(sp, Ω) do ω
        if ω.is_noise
            fix(omega, ω.value; force = true)
        else
            unfix(omega)
            fix(y[1].out, ω.value; force = true)
        end
    end
end
simulations = SDDP.simulate(model, 1, [:y]);
y = [stage[:y][1].out for stage in simulations[1]]
omega = NoiseTerm.(false, [1.0, 1.2, 1.4])
sampling_scheme = SDDP.Historical(tuple.(1:3, omega));
simulations = SDDP.simulate(model, 1, [:y]; sampling_scheme);
y = [stage[:y][1].out for stage in simulations[1]]
```
