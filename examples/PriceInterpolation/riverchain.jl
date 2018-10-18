#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    The code in this file runs the examples from the paper

    Downward, A., Dowson, O., and Baucke, R. (2018). On the convergence of a
    cutting plane method for multistage stochastic programming problems with
    stagewise dependent price uncertainty. Optimization Online.
=#

using JuMP, SDDP, Clp, Base.Test

struct PriceTurbine
    flowknots::Vector{Float64}
    powerknots::Vector{Float64}
end

struct PriceReservoir
    min::Float64
    max::Float64
    initial::Float64
    turbine::PriceTurbine
    spill_cost::Float64
    inflows::Vector{Float64}
end

function priceprocess(USE_AR1, lipschitz_constant = 1e6)
    b_t = [61.261, 56.716, 59.159, 66.080, 72.131, 76.708, 76.665, 76.071, 76.832, 69.970, 69.132, 67.176]
    alpha = 0.5
    beta = -0.5
    minprice = 40.0
    maxprice = 100.0
    pbar = 61.261
    noise = [0.5,1.5,2.5,3.5,4.5]
    NOISES = DiscreteDistribution(vcat(-reverse(noise), noise))

    function ar1price(price, noise, stage, markovstate)
        if stage > 1
            return alpha * price + (1-alpha) * b_t[stage] + noise
        else
            return price
        end
    end

    function ar2price(price, noise, stage, markovstate)
        # price is a Tuple{Float64, Float64}
        # price[1] is t-1, price[2] is t-2
        if stage > 1
            return (alpha * price[1] + (1-alpha) * b_t[stage] + beta * (price[1] - price[2]) + noise, price[1])
        else
            return price
        end
    end

    if USE_AR1
        return (t,i) -> DynamicPriceInterpolation(
            dynamics       = (p,w) -> ar1price(p, w, t, i),
            initial_price  = pbar,
            min_price      = minprice,
            max_price      = maxprice,
            noise          = NOISES,
            lipschitz_constant = lipschitz_constant
        )
    else
        return (t,i) -> DynamicPriceInterpolation(
            dynamics       = (p,w) -> ar2price(p, w, t, i),
            initial_price  = (pbar, pbar),
            min_price      = (minprice, minprice),
            max_price      = (maxprice, maxprice),
            noise          = NOISES,
            lipschitz_constant = lipschitz_constant
        )
    end
end

function buildmodel(USE_AR1, valley_chain, lipschitz_constant = 1e6)
    return SDDPModel(
                sense           = :Min,
                stages          = 12,
                objective_bound = -50_000.0,
                solver          = ClpSolver(),
                value_function   = priceprocess(USE_AR1, lipschitz_constant)
                                        ) do sp, stage
        N = length(valley_chain)
        turbine(i) = valley_chain[i].turbine
        @state(sp, valley_chain[r].min <= reservoir[r=1:N] <= valley_chain[r].max, reservoir0==valley_chain[r].initial)
        @variables(sp, begin
            70 >= outflow[r=1:N]      >= 0
            spill[r=1:N]        >= 0
            generation_quantity >= 0 # Total quantity of water
            # Proportion of levels to dispatch on
            0 <= dispatch[r=1:N, level=1:length(turbine(r).flowknots)] <= 1
        end)
        @constraints(sp, begin
            # flow from upper reservoir
            reservoir[1] == reservoir0[1] - outflow[1] - spill[1]
            # other flows
            flow[i=2:N], reservoir[i] == reservoir0[i] - outflow[i] - spill[i] + outflow[i-1] + spill[i-1]
            # Total quantity generated
            generation_quantity == sum(turbine(r).powerknots[level] * dispatch[r,level] for r in 1:N for level in 1:length(turbine(r).powerknots))
            # Flow out
            turbineflow[r=1:N], outflow[r] == sum(turbine(r).flowknots[level] * dispatch[r, level] for level in 1:length(turbine(r).flowknots))
            # Dispatch combination of levels
            dispatched[r=1:N], sum(dispatch[r, level] for level in 1:length(turbine(r).flowknots)) <= 1

            maxflow[i=1:N], outflow[i] <= reservoir0[i]
        end)
        if USE_AR1
            @stageobjective(sp,
                price -> sum(
                    valley_chain[i].spill_cost * spill[i] for i in 1:N
                    ) - price * generation_quantity
            )
        else
            @stageobjective(sp,
                price -> sum(
                    valley_chain[i].spill_cost * spill[i] for i in 1:N
                    ) - price[1] * generation_quantity
            )
        end
    end
end

function lipschitz_constant_experiment()
    river_chain = [
        # PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
        # PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
        # PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
        PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
        PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0])
    ]
    for alpha in [0.0, 10.0, 100.0, 1_000.0, 10_000.0, 100_000.0, 1_000_000.0]
        model = buildmodel(true, river_chain, alpha)
        solve(model,
            iteration_limit = 500,
            log_file="alpha_$(alpha).log")
    end
end

function runpaper(USE_AR1, valley_chain, name)
    m = buildmodel(USE_AR1, valley_chain)
    solve(m, iteration_limit = 2_000, log_file="$(name).log")

    SIMN = 1_000
    sim = simulate(m, SIMN, [:reservoir, :price, :generation_quantity])
    plt = SDDP.newplot()
    if USE_AR1
        SDDP.addplot!(plt, 1:SIMN, 1:12, (i, t)->sim[i][:price][t], title="Simulated Price", ylabel="Price (\$/MWH)")
    else
        SDDP.addplot!(plt, 1:SIMN, 1:12, (i, t)->sim[i][:price][t][1], title="Simulated Price", ylabel="Price (\$/MWH)")
    end
    SDDP.addplot!(plt, 1:SIMN, 1:12, (i, t)->sim[i][:generation_quantity][t], title="Offer", ylabel="Quantity (MWH)")
    for r in 1:length(valley_chain)
        SDDP.addplot!(plt, 1:SIMN, 1:12, (i, t)->sim[i][:reservoir][t][r], title="Storage (Reservoir $(r))", ylabel="Volume (m^3)")
    end
    SDDP.show("$(name).html", plt)
end

#=
    These were the commands used to run the examples in the paper
=#
# runpaper(
#     true,
#     [
#         PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
#         PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0])
#     ],
#     "example_one"
# )

# runpaper(
#     false,
#     [
#         PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
#         PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
#         PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
#         PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
#         PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0])
#     ],
#     "example_two"
# )

srand(123)
m = buildmodel(false, [
        PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
        PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0])
    ]
)
srand(123)
solve(m, iteration_limit = 20, print_level=0, cut_output_file="river.cuts")
@test getbound(m) >= -40_000

# test cut round trip with multiple price dimensions
srand(123)
m2 = buildmodel(false, [
        PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0]),
        PriceReservoir(0, 200, 100, PriceTurbine([50, 60, 70], [55, 65, 70]), 1000, [0])
    ]
)
try
    SDDP.loadcuts!(m2, "river.cuts")
finally
    rm("river.cuts")
end
srand(123)
SDDP.solve(m, iteration_limit=1, print_level=0)
srand(123)
SDDP.solve(m2, iteration_limit=1, print_level=0)
@test isapprox(getbound(m2), getbound(m), atol=1e-4)
