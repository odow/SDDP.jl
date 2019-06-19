import GLPK, JuMP, SDDP
using JuMP: @variable, @constraint, @objective

# Sample randomly 5 scenarios
nscen = 5;
probs = ones(nscen)/nscen;
noises = randn(nscen);

pb = SDDP.PolicyGraph(SDDP.LinearGraph(5),
            bellman_function = SDDP.BellmanFunction(lower_bound = 0.0),
            optimizer = JuMP.with_optimizer(GLPK.Optimizer)
           ) do sp, t

    @SDDP.variable(sp, x, SDDP.State, initial_value = 2.0)

    # noise
    @variable(sp,ξ)
    # control
    @variable(sp,c,Bin)
    # aux
    @variable(sp,y)

    # dynamics
    @constraint(sp, x.out == x.in + ξ + 2*c - 1)
    # aux
    @constraint(sp, y >=  x.out)
    @constraint(sp, y >= -x.out)

    # Noise
    SDDP.parameterize(sp, noises, probs) do ω
        JuMP.fix(ξ, ω)
    end

    # objective
    @SDDP.stageobjective(sp, y)
end

SDDP.train(pb; iteration_limit=20)
