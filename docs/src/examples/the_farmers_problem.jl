# # The farmer's problem
#
# _This problem is taken from Section 1.1 of the book Birge, J. R., & Louveaux,
# F. (2011). Introduction to Stochastic Programming. New York, NY: Springer New
# York. The wording in the problem description section is taken verbatim._

# ## Problem description
#
# "Consider a European farmer who specializes in raising wheat, corn, and sugar
# beets on his 500 acres of land. During the winter, [they want] to decide how
# much land to devote to each crop."

CROPS = [:wheat, :corn, :sugar_beet]

# "The farmer knows that at least 200 tons (T) of wheat and 240 T of corn are
# needed for cattle feed. These amounts can be raised on the farm or bought from
# a wholesaler. Any production in excess of the feeding requirement would be
# sold."

MIN_QUANTITIES = Dict(
    :wheat => 200.0,
    :corn => 240.0
)

# "Over the last decade, mean selling prices have been $170 and $150 per ton of
# wheat and corn, respectively. The purchase prices are 40% more than this due
# to the wholesaler’s margin and transportation costs."

SELL_PRICE = Dict(
    :wheat => 170.0,
    :corn => 150.0
)

BUY_PRICE = Dict(key => 1.4 * value for (key, value) in SELL_PRICE)

# "Another profitable crop is sugar beet, which [they expect] to sell at $36/T;
# however, the European Commission imposes a quota on sugar beet production. Any
# amount in excess of the quota can be sold only at $10/T. The farmer’s quota
# for next year is 6000 T."

SELL_PRICE[:sugar_beet] = 36.0

SELL_NO_QUOTA = copy(SELL_PRICE)
SELL_NO_QUOTA[:sugar_beet] = 10.0

QUOTA_MAX = Dict(
    :sugar_beet => 6_000.0
)

# "Based on past experience, the farmer knows that the mean yield on [their]
# land is roughly 2.5 T, 3 T, and 20 T per acre for wheat, corn, and sugar
# beets, respectively."

MEAN_YIELD = Dict(
    :wheat => 2.5,
    :corn => 3.0,
    :sugar_beet => 20.0
)

# To introduce uncertainty, "assume some correlation among the yields of the
# different crops. A very simplified representation of this would be to assume
# that years are good, fair, or bad for all crops, resulting in above average,
# average, or below average yields for all crops. To fix these ideas, _above_
# and _below_ average indicate a yield 20% above or below the mean yield."

YIELD_MULTIPLIER = Dict(
    :good => 1.2,
    :fair => 1.0,
    :bad => 0.8
)

# ## Mathematical description

# ## SDDP.jl code

# First up, load `SDDP.jl` and a solver. For this example, we use [`GLPK.jl`](https://github.com/JuliaOpt/GLPK.jl).

using SDDP, GLPK

# ### Common elements

function add_state_variables(subproblem)
    @variable(subproblem, area[c = CROPS] >= 0, SDDP.State, initial_value=0)
end

# ### First stage problem

function create_first_stage_problem(subproblem)
    area = subproblem[:area]
    @constraint(subproblem, sum(area[c].out for c in CROPS) <= 500.0)
    @stageobjective(subproblem, 0.0)
end

# ### Second stage problem

function second_stage_variables(subproblem)
    @variables(subproblem, begin
        yield[c=CROPS]
        buy[c=CROPS] >= 0
        sell_in_quota[c=CROPS] >= 0
        sell_no_quota[c=CROPS] >= 0
    end)
end

# #### Constraints

function second_stage_constraint_min_quantity(subproblem)
    @constraint(subproblem, [c=CROPS],
        yield[c].out + buy[c] - sell_in_quota[c] - sell_no_quota[c] >=
            get(MIN_QUANTITIES, c, 0.0)
    )
end

#

function second_stage_constraint_quota(subproblem)
    @constraint(subproblem, [c=CROPS],
        subproblem[:sell_in_quota][c] <= get(QUOTA_MAX, c, Inf)
    )
end

#

function second_stage_objective(subproblem)
    @stageobjective(subproblem,
        sum(
            SELL_PRICE[c] * subproblem[:sell_in_quota][c] +
            SELL_NO_QUOTA[c] * subproblem[:sell_no_quota][c] -
            BUY_PRICE[c] * subproblem[:buy][c]
        for c in CROPS)
    )
end

#

function second_stage_uncertainty(subproblem)
    SDDP.parameterize(subproblem, [:good, :fair, :bad]) do ω
        yield_multiplier = YIELD_MULTIPLIER[ω]
        for c in CROPS
            JuMP.fix(
                subproblem[:yield][c],
                MEAN_YIELD[c] * yield_multiplier * subproblem[:area][c].in
            )
        end
    end
end

#

function create_second_stage_problem(subproblem)
    second_stage_variables(subproblem)
    second_stage_constraint_quota(subproblem)
    second_stage_constraint_min_quantity(subproblem)
    second_stage_uncertainty(subproblem)
    second_stage_objective(subproblem)
end

# ### Putting it all together

model = SDDP.LinearPolicyGraph(
            stages = 2, lower_bound = 0.0,
            optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, stage
    add_state_variables(subproblem)
    if stage == 1
        create_first_stage_problem(subproblem)
    else
        create_second_stage_problem(subproblem)
    end
end

# ## Training a policy

SDDP.train(model; iteration_limit = 10)
