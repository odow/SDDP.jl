#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

# Add all available procs for parallel speeeeeeed
# addprocs(Sys.CPU_CORES-1)

# Load our favourite packages
using SDDP, JuMP, Gurobi#, CPLEX

# These commands get run on all processors
@everywhere begin

    # Intialise some random seed.
    # Note: You really shouldn't have the same seed on all processors or
    #   you end up doing the same thing N times!
    srand(1000*myid())
    # srand(1000)

    include("posm_data.jl")

    # Extract scenario data
    YEARS = unique(NIWA[:, 1])
    function EP(year, week)
        for i in 1:size(NIWA, 1)
            if NIWA[i, 1] == year && NIWA[i, 2] == week
                return NIWA[i, 3]
            end
        end
        return 0.0
    end
    function Rainfall(year, week)
        for i in 1:size(NIWA, 1)
            if NIWA[i, 1] == year && NIWA[i, 2] == week
                return NIWA[i, 4]
            end
        end
        return 0.0
    end

    # A helper set
    feed_type = [:supplement, :pasture]

    STATES = Dict(
        :SoilMoisture  => Dict(
            :min =>0,
            :max=>150,
            :initial=>150
        ),
        :PastureCover  => Dict(
            :min =>0,
            :max=>3500,
            :initial=>2200,
            :endminimum=>2200
        ),
        :FeedReserves  => Dict(
            :min =>0,
            :max=>50,
            :initial=>0,
            :finalvalue=>400
        ),
        :NumberMilking => Dict(
            :min =>0,
            :max=>175,
            :initial=>175
        )
    )

    # Helper functions for grass growth
    #   y = a * x * (1-x/N)
    function grwth(t, g)
        ep_max = maximum([EP(year, t) for year in YEARS])
        delta_max = 7 * (MOIR_INTERCEPT + MOIR_SLOPE * ep_max)
        4 * delta_max / STATES[:PastureCover][:max] * g * (1 - g / STATES[:PastureCover][:max])
    end
    #   dy/dx = a - 2ax/N
    function dGrwthdG(t, g)
        ep_max = maximum([EP(year, t) for year in YEARS])
        delta_max = 7 * (MOIR_INTERCEPT + MOIR_SLOPE * ep_max)
        4 * delta_max / STATES[:PastureCover][:max]  - 8 * delta_max / STATES[:PastureCover][:max]^2 * g
    end

end

m = SDDPModel(
    sense             = :Max,
    stages            = T,
    solver            = GurobiSolver(OutputFlag=0), # Note: I get some weird issues with Clp sometimes
    # solver            = CplexSolver(CPX_PARAM_SCRIND=0),
    objective_bound   = MAX_PROFIT
                                                        ) do sp, t, markovstate
    #--------------------------------------------------------------------------
    #   State Variables
    statekeys = collect(keys(STATES))
    @state(sp, STATES[s][:min] <= S[s=statekeys] <= STATES[s][:max], S0==STATES[s][:initial])

    @variables(sp, begin
        #----------------------------------------------------------------------
        #   Actions
        dry_off            >= 0 # number cows to dry off
        0 <= harvest <= MAX_HARVEST # tonnes
        buy_feed           >= 0 # tonnes
        feed[feed_type]    >= 0 # tonnes

        # Introduce a penalty variable to penalise
        #   when more than 8kg/day supplement is fed
        supplement_penalty >= 0  # kg

        evapotranspiration >= 0 # Actual Evapotranspiration (mm/week)
        drainage           >= 0 # Drainage slack variable (mm/week)
        pasture_growth     >= 0 # Growth (kgDM/Ha/week)
        rainfall           >= 0 # Rainfall (mm/week)
    end)

    #--------------------------------------------------------------------------
    # Stochastic constraints
    @scenarios(sp, year=YEARS, begin
        #   Rainfall
        rainfall <= 7*Rainfall(year, t)
        #   Potential evapotranspiration
        evapotranspiration <= 7 * EP(year, t)
    end)

    @constraints(sp, begin
        # Water limitations
        evapotranspiration <= S0[:SoilMoisture] + rainfall
        # Moir Function
        pasture_growth <= 7 * (MOIR_INTERCEPT + MOIR_SLOPE * evapotranspiration/7)

        #--------------------------------------------------------------------------
        #   Bounds
        # cannot dry more cows than milking
        dry_off        <= S0[:NumberMilking]
        # max growth
        pasture_growth <= (STATES[:PastureCover][:max] - S0[:PastureCover])

        #--------------------------------------------------------------------------
        #   State Conservation
        # Milking cow conservation
        S[:NumberMilking] <= S0[:NumberMilking] - dry_off

        # Pasture Conservation
        S[:PastureCover] <= S0[:PastureCover] +
            pasture_growth -             #   Pasture Growth
            1000 / AREA * ( #   Tonnes -> kg/Ha conversion
                feed[:pasture] +         #   Pasture fed to cows
                harvest                  #   Pasture turned into hay
            )

        # Supplement conservation
        S[:FeedReserves] <= S0[:FeedReserves] +
            buy_feed +         #   Quantity bought
            harvest -          #   Quantity harvested from pasture
            feed[:supplement]  #   Quantity fed out

        # Water conservation
        S[:SoilMoisture] == S0[:SoilMoisture] +
            rainfall - drainage - evapotranspiration
        # Drainage
        drainage >= (S0[:SoilMoisture] + rainfall - evapotranspiration) - 150

        #--------------------------------------------------------------------------
        #   Feed Conservation
        # Limited quantity of pasture
        feed[:pasture] <= S0[:PastureCover] * AREA / 1000. - harvest

        # Limited quantity of supplement
        feed[:supplement] <= S0[:FeedReserves] + buy_feed + harvest

        # Maximum supplement per day
        feed[:supplement] <= 8 * 7 * STATES[:NumberMilking][:max] / 1000.

        # Balance feed budget
        1000 * sum(feed[ft] for ft in feed_type) + supplement_penalty >=
            FEED[t, 1] * 7 * S[:NumberMilking] +
            FEED[t, 2] * 7 * (STATES[:NumberMilking][:max] - S[:NumberMilking]) +
            (BCS[t+1] - BCS[t]) * BCS_RESPONSE * STATES[:NumberMilking][:max]
    end)

    @constraint(sp, populationgrowth[g=linspace(0, STATES[:PastureCover][:max], 20)],
        pasture_growth <= grwth(t, g) + dGrwthdG(t, g) * (S0[:PastureCover]-g)
    )

    # Dry off constraint
    t > MAX_WEEK && @constraint(sp, S[:NumberMilking] == 0)

    #--------------------------------------------------------------------------
    #   Profit definition
    supplement_penalty_cost = 1 # $/kg above 8kg/day
    weekly_profit = 6.0 * MILK_SOLIDS[t] * 7 * S[:NumberMilking] -
        SUPPLEMENT_COST[t] * buy_feed -
        HARVEST_COST * harvest -
        1e-4 * drainage - # Note: this is to ensure we maintain water
        supplement_penalty_cost * supplement_penalty

    if t == T
        @variable(sp,  end_pasture_slack>=0)
        @constraint(sp, S[:PastureCover] + end_pasture_slack >= STATES[:PastureCover][:endminimum])
        append!(weekly_profit, S[:FeedReserves] * STATES[:FeedReserves][:finalvalue] - 1e6 * end_pasture_slack)
    end

    stageobjective!(sp, weekly_profit)

end
#
@time status = solve(m,
    max_iterations=150,
    simulation     = MonteCarloSimulation(
                        frequency = 10,
                        min       = 500,
                        max       = 500,
                        step      = 100
                             ),
    solve_type = Asyncronous()
)
#
# @time results = simulate(m, 5000, [:S, :feed], parallel=true)
#
# x = results[:Objective]
# Objective["Max"] = maximum(x)
# Objective["UQ"] = quantile(x, 0.75)
# Objective["Med"] = quantile(x, 0.5)
# Objective["LQ"] = quantile(x, 0.25)
# Objective["Min"] = minimum(x)
#
# for t=1:52
#     for (ql, qt) in [("Q1", 0.05), ("M", 0.5), ("Q2", 0.95)]
#         for state in STATES
#             Solution[t, state, ql] = quantile([s[state] for s in results[:S][t]], qt)
#         end
#         Solution[t, "pasture", ql] = quantile([r[:pasture] for r in results[:feed][t]], qt)
#         Solution[t, "supplement", ql] = quantile([r[:supplement] for r in results[:feed][t]], qt)
#     end
# end
