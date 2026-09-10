#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Here-and-now and hazard-decision

# SDDP.jl assumes that the agent gets to make a decision _after_ observing the
# realization of the random variable. This is called a _wait-and-see_ or
# _hazard-decision_ model. In contrast, you might want your agent to make
# decisions _before_ observing the random variable. This is called a
# _here-and-now_ or _decision-hazard_ model.

# !!! info
#     The terms decision-hazard and hazard-decision from the French _hasard_,
#     meaning chance. It could also have been translated as uncertainty-decision
#     and decision-uncertainty, but the community seems to have settled on the
#     transliteration _hazard_ instead. We like the hazard-decision and
#     decision-hazard terms because they clearly communicate the order of the
#     decision and the uncertainty.

# The purpose of this tutorial is to demonstrate how to model here-and-now
# decisions in SDDP.jl.

# This tutorial uses the following packages:

using SDDP
import HiGHS

# ## Hazard-decision formulation

# As an example, we're going to build a standard hydro-thermal scheduling
# model, with a single hydro-reservoir and a single thermal generation plant.
# In each of the four stages, we need to choose some mix of `u_thermal` and
# `u_hydro` to meet a demand of `9` units, where unmet demand is penalized at a
# rate of \$500/unit.

hazard_decision = SDDP.LinearPolicyGraph(;
    stages = 4,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, node
    @variables(sp, begin
        0 <= x_storage <= 8, (SDDP.State, initial_value = 6)
        u_thermal >= 0
        u_hydro >= 0
        u_unmet_demand >= 0
    end)
    @constraint(sp, u_thermal + u_hydro == 9 - u_unmet_demand)
    @constraint(sp, c_balance, x_storage.out == x_storage.in - u_hydro + 0)
    SDDP.parameterize(sp, [2, 3]) do ω_inflow
        return set_normalized_rhs(c_balance, ω_inflow)
    end
    @stageobjective(sp, 500 * u_unmet_demand + 20 * u_thermal)
end

# ## Decision-hazard formulation

# In the wait-and-see formulation, we get to decide the generation variables
# _after_ observing the realization of `ω_inflow`. However, a common modeling
# situation is that we need to decide the level of thermal generation
# `u_thermal` _before_ observing the inflow.

# SDDP.jl can model here-and-now decisions with a modeling trick: a wait-and-see
# decision in stage _t-1_ is equivalent to a here-and-now decision in stage _t_.

# In other words, we need to convert the `u_thermal` decision from a control
# variable that is decided in stage `t`, to a state variable that is decided in
# stage `t-1`. Here's our new model, with the three lines that have changed:

decision_hazard = SDDP.LinearPolicyGraph(;
    stages = 4,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, node
    @variables(sp, begin
        0 <= x_storage <= 8, (SDDP.State, initial_value = 6)
        u_thermal >= 0, (SDDP.State, initial_value = 0)  # <-- changed
        u_hydro >= 0
        u_unmet_demand >= 0
    end)
    @constraint(sp, u_thermal.in + u_hydro == 9 - u_unmet_demand)  # <-- changed
    @constraint(sp, c_balance, x_storage.out == x_storage.in - u_hydro + 0)
    SDDP.parameterize(sp, [2, 3]) do ω
        return set_normalized_rhs(c_balance, ω)
    end
    @stageobjective(sp, 500 * u_unmet_demand + 20 * u_thermal.in) # <-- changed
end

# Can you understand the reformulation? In each stage, we now use the value of
# `u_thermal.in` instead of `u_thermal`, and the value of the outgoing
# `u_thermal.out` is the here-and-how decision for stage `t+1`.

# (If you can spot a "mistake" with this model, don't worry, we'll fix it below.
# Presenting it like this simplifies the exposition.)

# ## Comparison

# Let's compare the cost of operating the two models:

function train_and_compute_cost(model)
    SDDP.train(model; print_level = 0)
    return println("Cost = \$", SDDP.calculate_bound(model))
end

train_and_compute_cost(hazard_decision)

#-

train_and_compute_cost(decision_hazard)

# This suggests that choosing the thermal generation before observing the inflow
# adds a cost of \$250. But does this make sense?

# If we look carefully at our `decision_hazard` model, the incoming value of
# `u_thermal.in` in the first stage is fixed to the `initial_value` of `0`.
# Therefore, we must always meet the full demand with `u_hydro`, which we cannot
# do without incurring unmet demand.

# To allow the model to choose an optimal level of `u_thermal` in the
# first-stage, we need to add an extra stage that is deterministic with no
# stage objective.

# ## Fixing the decision-hazard

# In the following model, we now have five stages, so that stage `t+1` in
# `decision_hazard_2` corresponds to stage `t` in `decision_hazard`. We've also
# added an `if`-statement, which adds different constraints depending on the
# node. Note that we need to add an `x_storage.out == x_storage.in` constraint
# because the storage can't change in this new first-stage.

decision_hazard_2 = SDDP.LinearPolicyGraph(;
    stages = 5,  # <-- changed
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, node
    @variables(sp, begin
        0 <= x_storage <= 8, (SDDP.State, initial_value = 6)
        u_thermal >= 0, (SDDP.State, initial_value = 0)
        u_hydro >= 0
        u_unmet_demand >= 0
    end)
    if node == 1                                        # <-- new
        @constraint(sp, x_storage.out == x_storage.in)  # <-- new
        @stageobjective(sp, 0)                          # <-- new
    else
        @constraint(sp, u_thermal.in + u_hydro == 9 - u_unmet_demand)
        @constraint(sp, c_balance, x_storage.out == x_storage.in - u_hydro + 0)
        SDDP.parameterize(sp, [2, 3]) do ω
            return set_normalized_rhs(c_balance, ω)
        end
        @stageobjective(sp, 500 * u_unmet_demand + 20 * u_thermal.in)
    end
end

train_and_compute_cost(decision_hazard_2)

# Now we find that the cost of choosing the thermal generation before observing
# the inflow adds a much more reasonable cost of \$10.

# ## Summary

# To summarize, the difference between here-and-now and wait-and-see variables
# is a modeling choice.

# To create a here-and-now decision, add it as a state variable to the
# previous stage

# In some cases, you'll need to add an additional "first-stage" problem to
# enable the model to choose an optimal value for the here-and-now decision
# variable. You do not need to do this if the first stage is deterministic. You
# must make sure that the subproblem is feasible for all possible incoming
# values of the here-and-now decision variable.
