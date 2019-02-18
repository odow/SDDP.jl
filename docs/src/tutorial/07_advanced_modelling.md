# Basic VII: modelling tips

This tutorial discusses some different modelling tips.

```@meta
DocTestSetup = quote
    using Kokako, GLPK
end
```

## Multi-dimensional state variables

Just like normal JuMP variables, it is possible to create containers of state
variables.

```jldoctest; filter=r"A policy graph.+"s
Kokako.LinearPolicyGraph(
        stages=1, lower_bound = 0, optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, t
    # A scalar state variable.
    @variable(subproblem, x >= 0, Kokako.State, initial_value = 0)
    println("Lower bound of outgoing x is: ", JuMP.lower_bound(x.out))

    # A vector of state variables.
    @variable(subproblem, y[i = 1:2] >= i, Kokako.State, initial_value = i)
    println("Lower bound of outgoing y[1] is: ", JuMP.lower_bound(y[1].out))

    # A JuMP.Containers.DenseAxisArray of state variables.
    @variable(subproblem,
        z[i = 3:4, j = [:A, :B]] >= i, Kokako.State, initial_value = i)
    println("Lower bound of outgoing z[3, :B] is: ", JuMP.lower_bound(z[3, :B].out))
end

nothing # hide

# output

Lower bound of outgoing x is: 0.0
Lower bound of outgoing y[1] is: 1.0
Lower bound of outgoing z[3, :B] is: 3.0
```

## Multi-dimensional noise terms

Multi-dimensional stagewise-independent random variables can be created by
forming the Cartesian product of the random variables.

```jldoctest; filter=[r"\(value = \d, coefficient = \d\)", r"1\-element.+"s]
model = Kokako.LinearPolicyGraph(
        stages=3, lower_bound = 0, optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, t
    @variable(subproblem, x, Kokako.State, initial_value = 0.0)

    support = [(value = v, coefficient = c) for v in [1, 2] for c in [3, 4, 5]]
    probability = [pv * pc for pv in [0.5, 0.5] for pc in [0.3, 0.5, 0.2]]
    Kokako.parameterize(subproblem, support, probability) do ω
        JuMP.fix(x.out, ω.value)
        @stageobjective(subproblem, ω.coefficient * x.out)
        println("ω is: ", ω)
    end
end

Kokako.simulate(model, 1)

nothing # hide

# output

ω is: (value = 1, coefficient = 4)
ω is: (value = 1, coefficient = 3)
ω is: (value = 2, coefficient = 4)
```

This concludes or series of basic introductory tutorials for `Kokako.jl`. When
you're ready, continue to our intermediate series of tutorials, beginning with
[Intermediate I: risk](@ref).
