#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Duality handlers

# The purpose of this tutorial is to demonstrate trivial examples that expose
# the strengths and weaknesses of each duality handler.

# For more information on SDDP.jl's duality handlers, see [Integrality](@ref).

# This tutorial uses the following packages:

using SDDP
import HiGHS

# Note that these trivial examples exposed a bug in HiGHS, which we work-around
# by turning off presolve. In real examples, you should probably leave the
# default options.

Optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")

# ## Training function

# First, we need a function to simplify our testing:

function train_and_evaluate_bounds(
    model_fn::Function,
    duality_handler::SDDP.AbstractDualityHandler;
    print_level::Int = 0,
    kwargs...,
)
    model = model_fn()
    SDDP.train(model; print_level, duality_handler, kwargs...)
    simulations = SDDP.simulate(model, 1)
    lower_bound = SDDP.calculate_bound(model)
    println("lower_bound: $(lower_bound)")
    upper_bound = sum(data[:stage_objective] for data in only(simulations))
    println("upper_bound: $(upper_bound)")
    return
end

# This function builds a new model, trains it using the provided `duality_handler`,
# and then prints the lower and upper bounds. We'll assume that the models we
# pass in are deterministic so we need to conduct only a single simulation to
# evaluate the upper bound.

# !!! danger
#     The most important thing to keep in mind when reading this tutorial is
#     that SDDP.jl is not guaranteed to find a globally optimal policy. No
#     matter what options we select, there may be a gap between the lower and
#     upper bound.

# ## [`ContinuousConicDuality`](@id section_continuous)

# The default duality handler in SDDP.jl is [`ContinuousConicDuality`](@ref).
# To compute a cut, it solves the continuous relaxation of the MIP.

# In the same way that solution to a relaxed linear program may be far from the
# optimal MIP solution, the biggest downside to [`ContinuousConicDuality`](@ref)
# is that many models have large gaps between the lower and upper bound:

function model_1()
    return SDDP.LinearPolicyGraph(;
        stages = 2,
        lower_bound = 0.0,
        optimizer = Optimizer,
    ) do sp, t
        @variable(sp, x, Bin, SDDP.State, initial_value = 1.0)
        @variable(sp, y, Bin)
        @constraint(sp, x.out == x.in)
        if t == 1
            @stageobjective(sp, x.out)
        else
            @stageobjective(sp, y)
            @constraint(sp, y >= x.in - 0.5)
        end
    end
end

train_and_evaluate_bounds(model_1, SDDP.ContinuousConicDuality())

# ## [`StrengthenedConicDuality`](@id section_strengthened)

# One technique to improve upon [`ContinuousConicDuality`](@ref) is
# [`StrengthenedConicDuality`](@ref). Without going into the technical details,
# both use the continuous relaxation to compute a valid subgradient for the cut.
# [`StrengthenedConicDuality`](@ref) then tries to improve the cut by solving an
# additional integer program.

# In this example, [`StrengthenedConicDuality`](@ref) can improve upon the lower
# bound and prove that the upper bound of `2.0` is optimal:

train_and_evaluate_bounds(model_1, SDDP.StrengthenedConicDuality())

# Sometimes, however, [`StrengthenedConicDuality`](@ref) cannot improve upon
# [`ContinuousConicDuality`](@ref):

function model_2()
    return SDDP.LinearPolicyGraph(;
        stages = 2,
        lower_bound = 0.0,
        optimizer = Optimizer,
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 0.1)
        @variable(sp, y, Int)
        @constraint(sp, x.out == x.in)
        if t == 1
            @stageobjective(sp, x.out)
        else
            @stageobjective(sp, y)
            @constraint(sp, y >= x.in + 0.1)
            @constraint(sp, y >= -x.in + 0.1)
        end
    end
end

train_and_evaluate_bounds(model_2, SDDP.ContinuousConicDuality())

#-

train_and_evaluate_bounds(model_2, SDDP.StrengthenedConicDuality())

# Even though it is sometimes tighter than [`ContinuousConicDuality`](@ref) and
# it can never be worse, [`StrengthenedConicDuality`](@ref) is not the default
# duality handler because it is more expensive to compute; it solves a
# mixed-integer program whereas [`ContinuousConicDuality`](@ref) solves a
# continuous relaxation.

# ## [`LagrangianDuality`](@id section_lagrangian)

# A technique to improve upon [`StrengthenedConicDuality`](@ref) is
# [`LagrangianDuality`](@ref). Without going into the technical details,
# both use the continuous relaxation to compute a valid subgradient for the cut,
# but, where [`StrengthenedConicDuality`](@ref) tries to improve the cut by
# solving a single additional integer program, [`LagrangianDuality`](@ref) may
# solve many integer programs.

# [`LagrangianDuality`](@ref) finds the optimal policy for `model_1`:

train_and_evaluate_bounds(model_1, SDDP.LagrangianDuality())

# and also for `model_2`:

train_and_evaluate_bounds(model_2, SDDP.LagrangianDuality())

# Sometimes, however, [`LagrangianDuality`](@ref) does not close the gap:

function model_3()
    return SDDP.LinearPolicyGraph(;
        stages = 2,
        lower_bound = -1.0,
        optimizer = Optimizer,
    ) do sp, t
        @variable(sp, -1 <= x <= 0.5, SDDP.State, initial_value = 0.0)
        @variable(sp, y)
        @stageobjective(sp, y)
        if t == 1
            @constraint(sp, y >= x.out)
            @constraint(sp, y >= -x.out)
        else
            @variable(sp, z, Bin)
            @constraint(sp, y >= 1 - x.in - 3 * z)
            @constraint(sp, y >= 1 + x.in - 3 * (1 - z))
        end
    end
end

train_and_evaluate_bounds(model_3, SDDP.LagrangianDuality())

# but it may still be better than [`StrengthenedConicDuality`](@ref):

train_and_evaluate_bounds(model_3, SDDP.StrengthenedConicDuality())

# The algorithm behind [`LagrangianDuality`](@ref) is significantly more
# complicated than [`StrengthenedConicDuality`](@ref). For some models it can
# be helpful, for others, the increased computational cost is not worth the
# improvement in the tightness of the value function.

# ## Different policies

# So far, the different duality handlers have led to different lower bounds, but
# identical upper bounds. This is an artifact of our trivial examples. Using a
# more sophisticated duality handler can improve the lower bound _and_ lead to
# a cheaper policy (the upper bound). Here's an example:

function model_4()
    return SDDP.LinearPolicyGraph(;
        stages = 2,
        lower_bound = -1.0,
        optimizer = Optimizer,
    ) do sp, t
        @variable(sp, -1 <= x <= 0.5, SDDP.State, initial_value = 0.0)
        @variable(sp, y)
        @variable(sp, z, Bin)
        if t == 1
            @stageobjective(sp, -0.1 * x.out)
        else
            @stageobjective(sp, y)
            @constraint(sp, y >= 1 - x.in - 3 * z)
            @constraint(sp, y >= 1 + x.in - 3 * (1 - z))
        end
    end
end

# [`ContinuousConicDuality`](@ref) finds in a policy that costs `0.45`:

train_and_evaluate_bounds(model_4, SDDP.ContinuousConicDuality())

# whereas [`LagrangianDuality`](@ref) finds a policy that costs `0.1`:

train_and_evaluate_bounds(model_4, SDDP.LagrangianDuality())

# This relationship is not guaranteed to hold. In some models
# [`ContinuousConicDuality`](@ref) may find a cheaper policy than
# [`LagrangianDuality`](@ref), even though the latter finds a tighter lower
# bound. In general, you should experiment with different duality handlers to
# see what works best for your problem.

# ## [`BanditDuality`](@id section_bandit)

# The trade-off between the computational cost and the tightness of a
# formulation can be tricky to manage. SDDP.jl includes [`BanditDuality`](@ref),
# which is an algorithm that does not appear in the published academic
# literature. The [`BanditDuality`](@ref) duality handler treats the problem of
# choosing a duality handler for each iteration of the SDDP algorithm as a
# multi-armed bandit problem, where the reward is the change in the lower bound
# after each iteration per second of computation time. The multi-armed bandit
# problem allows us to trade off many fast but weak iterations of
# [`ContinuousConicDuality`](@ref) against a small number of relatively strong
# iterations of [`LagrangianDuality`](@ref).

duality_handler = SDDP.BanditDuality(
    SDDP.ContinuousConicDuality(),
    SDDP.StrengthenedConicDuality(),
    SDDP.LagrangianDuality(),
)
train_and_evaluate_bounds(model_1, duality_handler)

#-

train_and_evaluate_bounds(model_2, duality_handler)

#-

train_and_evaluate_bounds(model_3, duality_handler)

# The [`BanditDuality`](@ref) is often a very good choice to use in practice.
# It is not the default because a tighter lower bound does not always lead to a
# better policy, so we opt for the simplest and fastest duality handler as the
# default.

# ## [`FixedDiscreteDuality`](@id section_fixed)

# An alternative to [`StrengthenedConicDuality`](@ref) is [`FixedDiscreteDuality`](@ref).

# It works by first solving the mixed-integer problem, fixing the discrete
# variables to their optimal value, and then solving the continuous relaxation.
# The cut from the continuous relaxation is then modified by solving another
# mixed-integer problem to ensure that it remains globally valid.

# For some models, [`FixedDiscreteDuality`](@ref) can find solutions that are
# tighter than [`StrengthenedConicDuality`](@ref):

train_and_evaluate_bounds(model_2, SDDP.FixedDiscreteDuality())

#-

train_and_evaluate_bounds(model_2, SDDP.StrengthenedConicDuality())

# Othertimes, it is weaker:

train_and_evaluate_bounds(model_1, SDDP.FixedDiscreteDuality())

#-

train_and_evaluate_bounds(model_1, SDDP.StrengthenedConicDuality())
