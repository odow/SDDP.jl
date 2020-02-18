# # The farmer's problem
#
# _This problem is taken from Section 1.1 of the book Birge, J. R., & Louveaux,
# F. (2011). Introduction to Stochastic Programming. New York, NY: Springer New
# York. Paragraphs in quotes are taken verbatim._

# ## Problem description

# > Consider a European farmer who specializes in raising wheat, corn, and sugar
# > beets on his 500 acres of land. During the winter, [they want] to decide how
# > much land to devote to each crop.

# > The farmer knows that at least 200 tons (T) of wheat and 240 T of corn are
# > needed for cattle feed. These amounts can be raised on the farm or bought
# > from a wholesaler. Any production in excess of the feeding requirement would
# > be sold.

# > Over the last decade, mean selling prices have been \\\$170 and \\\$150 per
# > ton of wheat and corn, respectively. The purchase prices are 40% more than
# > this due to the wholesaler’s margin and transportation costs.

# > Another profitable crop is sugar beet, which [they expect] to sell at
# > \\\$36/T; however, the European Commission imposes a quota on sugar beet
# > production. Any amount in excess of the quota can be sold only at \\\$10/T.
# > The farmer’s quota for next year is 6000 T."

# > Based on past experience, the farmer knows that the mean yield on [their]
# > land is roughly 2.5 T, 3 T, and 20 T per acre for wheat, corn, and sugar
# > beets, respectively.

# > [To introduce uncertainty,] assume some correlation among the yields of the
# > different crops. A very simplified representation of this would be to assume
# > that years are good, fair, or bad for all crops, resulting in above average,
# > average, or below average yields for all crops. To fix these ideas, _above_
# > and _below_ average indicate a yield 20% above or below the mean yield.

# ## Problem data

# The area of the farm.

MAX_AREA = 500.0

# There are three crops:

CROPS = [:wheat, :corn, :sugar_beet]

# Each of the crops has a different planting cost (\\\$/acre).

PLANTING_COST = Dict(
    :wheat      => 150.0,
    :corn       => 230.0,
    :sugar_beet => 260.0
)

# The farmer requires a minimum quantity of wheat and corn, but not of sugar
# beet (tonnes).

MIN_QUANTITIES = Dict(
    :wheat      => 200.0,
    :corn       => 240.0,
    :sugar_beet =>   0.0
)

# In Europe, there is a quota system for producing crops. The farmer owns the
# following quota for each crop (tonnes):

QUOTA_MAX = Dict(
    :wheat      =>     Inf,
    :corn       =>     Inf,
    :sugar_beet => 6_000.0
)

# The farmer can sell crops produced under the quota for the following amounts
# (\\\$/tonne):

SELL_IN_QUOTA = Dict(
    :wheat      => 170.0,
    :corn       => 150.0,
    :sugar_beet =>  36.0
)

# If they sell more than their alloted quota, the farmer earns the following on
# each tonne of crop above the quota (\\\$/tonne):

SELL_NO_QUOTA = Dict(
    :wheat      =>  0.0,
    :corn       =>  0.0,
    :sugar_beet => 10.0
)

# The purchase prices for wheat and corn are 40% more than their sales price.
# However, the description does not address the purchase price of sugar beet.
# Therefore, we use a large value of \\\$1,000/tonne.

BUY_PRICE = Dict(
    :wheat      =>   238.0,
    :corn       =>   210.0,
    :sugar_beet => 1_000.0
)

# On average, each crop has the following yield in tonnes/acre:

MEAN_YIELD = Dict(
    :wheat      =>  2.5,
    :corn       =>  3.0,
    :sugar_beet => 20.0
)

# However, the yield is random. In good years, the yield is +20% above average,
# and in bad years, the yield is -20% below average.

YIELD_MULTIPLIER = Dict(
    :good => 1.2,
    :fair => 1.0,
    :bad  => 0.8
)

# ## Mathematical formulation

# ## SDDP.jl code

# !!! note
#     In what follows, we make heavy use of the fact that you can look up
#     variables by their symbol name in a JuMP model as follows:
#     ```julia
#     @variable(model, x)
#     model[:x]
#     ```
#     Read the [JuMP documentation](http://www.juliaopt.org/JuMP.jl/v0.19/variables/)
#     if this isn't familiar to you.

# First up, load `SDDP.jl` and a solver. For this example, we use
# [`GLPK.jl`](https://github.com/JuliaOpt/GLPK.jl).

using SDDP, GLPK

# ### State variables

# State variables are the information that flows between stages. In our example,
# the state variables are the areas of land devoted to growing each crop.

function add_state_variables(subproblem)
    @variable(subproblem, area[c = CROPS] >= 0, SDDP.State, initial_value=0)
end

# ### First stage problem

# We can only plant a maximum of 500 acres, and we want to minimize the planting
# cost

function create_first_stage_problem(subproblem)
    @constraint(subproblem,
        sum(subproblem[:area][c].out for c in CROPS) <= MAX_AREA)
    @stageobjective(subproblem,
        -sum(PLANTING_COST[c] * subproblem[:area][c].out for c in CROPS))
end

# ### Second stage problem

# Now let's consider the second stage problem. This is more complicated than
# the first stage, so we've broken it down into four sections:
# 1) control variables
# 2) constraints
# 3) the objective
# 4) the uncertainty

# First, let's add the second stage control variables.

# #### Variables

# We add four types of control variables. Technically, the `yield` isn't a
# control variable. However, we add it as a dummy "helper" variable because it
# will be used when we add uncertainty.

function second_stage_variables(subproblem)
    @variables(subproblem, begin
        0 <= yield[c=CROPS]                          # tonnes/acre
        0 <= buy[c=CROPS]                            # tonnes
        0 <= sell_in_quota[c=CROPS] <= QUOTA_MAX[c]  # tonnes
        0 <= sell_no_quota[c=CROPS]                  # tonnes
    end)
end

# #### Constraints

# We need to define is the minimum quantity constraint. This ensures that
# `MIN_QUANTITIES[c]` of each crop is produced.

function second_stage_constraint_min_quantity(subproblem)
    @constraint(subproblem, [c=CROPS],
        subproblem[:yield][c] + subproblem[:buy][c] -
        subproblem[:sell_in_quota][c] - subproblem[:sell_no_quota][c] >=
        MIN_QUANTITIES[c])
end

# #### Objective

# The objective of the second stage is to maximise revenue from selling crops,
# less the cost of buying corn and wheat if necessary to meet the minimum
# quantity constraint.

function second_stage_objective(subproblem)
    @stageobjective(subproblem,
        sum(
            SELL_IN_QUOTA[c] * subproblem[:sell_in_quota][c] +
            SELL_NO_QUOTA[c] * subproblem[:sell_no_quota][c] -
            BUY_PRICE[c] * subproblem[:buy][c]
        for c in CROPS)
    )
end

# #### Random variables

# Then, in the [`SDDP.parameterize`](@ref) function, we set the coefficient
# using `JuMP.set_normalized_coefficient`.

function second_stage_uncertainty(subproblem)
    @constraint(subproblem, uncertainty[c=CROPS],
        1.0 * subproblem[:area][c].in == subproblem[:yield][c])
    SDDP.parameterize(subproblem, [:good, :fair, :bad]) do ω
        for c in CROPS
            JuMP.set_normalized_coefficient(
                uncertainty[c],
                subproblem[:area][c].in,
                MEAN_YIELD[c] * YIELD_MULTIPLIER[ω]
            )
        end
    end
end

# ### Putting it all together

# Now we're ready to build the multistage stochastic programming model. In
# addition to the things already discussed, we need a few extra pieces of
# information.
#
# First, we maximizing, so we set `sense = :Max`. Second, we need to provide a
# valid upper bound. (See [Choosing an initial bound](@ref) for more on this.)
# We know from Birge and Louveaux that the optimal solution is \\\$108,390.
# So, let's choose \\\$500,000 just to be safe.

# Here is the full model.

model = SDDP.LinearPolicyGraph(
            stages = 2,
            sense = :Max,
            upper_bound = 500_000.0,
            optimizer = GLPK.Optimizer,
            direct_mode = false
        ) do subproblem, stage
    add_state_variables(subproblem)
    if stage == 1
        create_first_stage_problem(subproblem)
    else
        second_stage_variables(subproblem)
        second_stage_constraint_min_quantity(subproblem)
        second_stage_uncertainty(subproblem)
        second_stage_objective(subproblem)
    end
end

# ## Training a policy

# Now that we've built a model, we need to train it using [`SDDP.train`](@ref).
# The keyword `iteration_limit` stops the training after 20 iterations. See
# [Choose a stopping rule](@ref) for other ways to stop the training.

SDDP.train(model; iteration_limit = 20)

# ## Checking the policy

# Birge and Louveaux report that the optimal objective value is \\\$108,390.
# Check that we got the correct solution using [`SDDP.calculate_bound`](@ref):

@assert SDDP.calculate_bound(model) == 108_390.0
