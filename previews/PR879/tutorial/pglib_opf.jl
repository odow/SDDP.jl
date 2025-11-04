#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Alternative forward models

# This example demonstrates how to train convex and non-convex models.

# This example uses the following packages:

using SDDP
import Ipopt
import PowerModels
import Test

# ## Formulation

# For our model, we build a simple optimal power flow model with a single
# hydro-electric generator.

# The formulation of our optimal power flow problem depends on `model_type`,
# which must be one of the `PowerModels` formulations.

# (To run locally, download [`pglib_opf_case5_pjm.m`](pglib_opf_case5_pjm.m) and
# update `filename` appropriately.)

function build_model(model_type)
    filename = joinpath(@__DIR__, "pglib_opf_case5_pjm.m")
    data = PowerModels.parse_file(filename)
    return SDDP.PolicyGraph(
        SDDP.UnicyclicGraph(0.95);
        sense = :Min,
        lower_bound = 0.0,
        optimizer = Ipopt.Optimizer,
    ) do sp, t
        power_model = PowerModels.instantiate_model(
            data,
            model_type,
            PowerModels.build_opf;
            jump_model = sp,
        )
        ## Now add hydro power models. Assume that generator 5 is hydro, and the
        ## rest are thermal.
        pg = power_model.var[:it][:pm][:nw][0][:pg][5]
        sp[:pg] = pg
        @variable(sp, x >= 0, SDDP.State, initial_value = 10.0)
        @variable(sp, deficit >= 0)
        @constraint(sp, balance, x.out == x.in - pg + deficit)
        @stageobjective(sp, objective_function(sp) + 1e6 * deficit)
        SDDP.parameterize(sp, [0, 2, 5]) do ω
            return SDDP.set_normalized_rhs(balance, ω)
        end
        return
    end
end

# ## Training a convex model

# We can build and train a convex approximation of the optimal power flow
# problem.

# The problem with the convex model is that it does not accurately simulate the
# true dynamics of the problem. Therefore, it under-estimates the true cost of
# operation.

convex = build_model(PowerModels.DCPPowerModel)
SDDP.train(convex; iteration_limit = 10)

# To more accurately simulate the dynamics of the problem, a common approach is
# to write the cuts representing the policy to a file, and then read them into
# a non-convex model:

SDDP.write_cuts_to_file(convex, "convex.cuts.json")
non_convex = build_model(PowerModels.ACPPowerModel)
SDDP.read_cuts_from_file(non_convex, "convex.cuts.json")

# Now we can simulate `non_convex` to evaluate the policy.

result = SDDP.simulate(non_convex, 1)

# A problem with reading and writing the cuts to file is that the cuts have been
# generated from trial points of the convex model. Therefore, the policy may be
# arbitrarily bad at points visited by the non-convex model.

# ## Training a non-convex model

# We can also build and train a non-convex formulation of the optimal power flow
# problem.

# The problem with the non-convex model is that because it is non-convex,
# SDDP.jl may find a sub-optimal policy. Therefore, it may over-estimate the
# true cost of operation.

non_convex = build_model(PowerModels.ACPPowerModel)
SDDP.train(non_convex; iteration_limit = 10)
result = SDDP.simulate(non_convex, 1)

# ## Combining convex and non-convex models

# To summarize, training with the convex model constructs cuts at points that
# may never be visited by the non-convex model, and training with the non-convex
# model may construct arbitrarily poor cuts because a key assumption of SDDP is
# convexity.

# As a compromise, we can train a policy using a combination of the convex and
# non-convex models; we'll use the non-convex model to generate trial points on
# the forward pass, and we'll use the convex model to build cuts on the backward
# pass.

convex = build_model(PowerModels.DCPPowerModel)

#-

non_convex = build_model(PowerModels.ACPPowerModel)

# To do so, we train `convex` using the [`SDDP.AlternativeForwardPass`](@ref)
# forward pass, which simulates the model using `non_convex`, and we use
# [`SDDP.AlternativePostIterationCallback`](@ref) as a post-iteration callback,
# which copies cuts from the `convex` model back into the `non_convex` model.

SDDP.train(
    convex;
    forward_pass = SDDP.AlternativeForwardPass(non_convex),
    post_iteration_callback = SDDP.AlternativePostIterationCallback(non_convex),
    iteration_limit = 10,
)

# In practice, if we were to simulate `non_convex` now, we should obtain a
# better policy than either of the two previous approaches.
