using GLPK
using JuMP
using SDDP
using Test

# easy case where there is a unique "tight" dual (current method won't guarantee we get it- a hack would be to set method to simplex)
model = SDDP.LinearPolicyGraph(
    stages = 2,
    sense = :Min,
    lower_bound = 0,
    optimizer = with_optimizer(GLPK.Optimizer)) do sp, stage

    @variable(sp, 0 <= x[i in 1:2] <= 1, SDDP.State, initial_value = 0)
    @variable(sp, y)

    if stage == 1
        @stageobjective(sp, 0)
    else
        @constraint(sp, y >= x[1].in + x[2].in)
        # @constraint(sp, x[1].in == 1) # TODO why didn't this work when JuMP.fix works?
        # @constraint(sp, x[2].in == 1)
        fix(x[1].in, 0)
        fix(x[2].in, 0)
        @stageobjective(sp, y)
    end
end
dual_vars = -2 * ones(2)
node = model.nodes[2]
obj = SDDP._kelley(node, dual_vars)
@test all(duals .>= -ones(2))
@test obj == 0.0


# 'exclusive or' function, no obvious choice of "tightest" dual
model = SDDP.LinearPolicyGraph(
    stages = 2,
    sense = :Min,
    lower_bound = 0,
    optimizer = with_optimizer(GLPK.Optimizer)) do sp, stage

    @variable(sp, 0 <= x[i in 1:2] <= 1, SDDP.State, initial_value = 1)
    @variable(sp, y)

    if stage == 1
        @stageobjective(sp, 0)
    else
        @constraints(sp, begin
            y >= x[1].in - x[2].in
            y >= x[2].in - x[1].in
            y <= x[1].in + x[2].in
            y <= 2 - x[1].in - x[2].in
        end)
        fix(x[1].in, 0)
        fix(x[2].in, 0)
        @stageobjective(sp, y)
    end
end
dual_vars = -2 * ones(2)
node = model.nodes[2]
obj = SDDP._kelley(node, dual_vars)
@test sum(duals) <= 1
@test obj == 0.0


;
