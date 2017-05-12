#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

#==

Example: The Widget Producer

A company produces perishable, divisible widgets (i.e. a continuous product,
rather than discrete units) for sale on a spot market each month. The quantity
of widgets they produce is uncertain and comes from some finite distribution.

The company can store the widgets, but, over time, the value of the widgets
decreases. Eventually the widgets perish and are worthless.

The spot price is determined by an auction system, and so varies from month to
month, but demonstates serial correlation. In each auction, there is sufficient
demand that the widget producer finds a buyer for all their widgets, regardless
of the quantity they supply. Furthermore, the spot price is independent of the
widget producer (they are a small player in the market).

The spot price is highly volatile, and is the result of a
process that is out of the control of the company. To counteract their price
risk, the company wishes to engage in a forward contracting programme.

This forward contracting programme is a deal for physical widgets at a future
date in time. The company can engage in contracts for sales two, three and six
months in the future.

The futures price is the current spot price, plus some forward contango (the
buyers gain certainty that they will receive the widgets in the future).

In general, the widget company should forward contract (since they reduce their
price risk), however they also have production risk. Therefore, it may be the
case that they forward contract a fixed amount, but find that they do not
produce enough widgets to meet the fixed demand. They are then forced to buy
additional widgets on the spot market.

The goal of the widget company is to choose the extent to which they forward
contract in order to maximise (risk-adjusted) revenues, whilst managing their
production risk.

==#

# These are a few of my favourite things...
using SDDP, JuMP, Gurobi

#==
    To enable the asyncronous solver, we need to copy the data we have to all
    the processors. To do that we can use the @everywhere macro.
==#
@everywhere begin

    #==
        For repeatability we should set the random number seed. However, if we
        choose the same seed on all processors, then we end up doing identical
        things! myid() returns an integer of the pocessor id, so this way we get
        a unique random seed for all processors.
    ==#
    srand(myid() * 10)

    #==
        The spot price is a log-normal AR(1) model:

            log(y(t)) = a + b * log(y(t-1)) + e,

        where a and b are constants, log(x) is the natural logarithm function, and
        e is a random noise.
    ==#
    const a = 1.01
    const b = 1.0
    const sampled_errors = [
        -0.1290, -0.1010, -0.0814, -0.0661, -0.0530, -0.0412,
        -0.0303, -0.0199, -0.00987, 0.0, 0.00987, 0.0199,
        0.0303, 0.0412, 0.0530, 0.0661, 0.0814, 0.1010, 0.1290]

    # A helper function to bind x in [lower, upper]
    box(x, lower, upper) = min(upper, max(lower, x))

    # The log AR(1) model
    price_dynamics(p, w, t, i) = box(a * exp(log(p) + w), 3, 9)

    # The spot price in period 0
    const initial_price = 4.50

    # There is a transaction cost to trading in the forward contracts
    const transaction_cost = 0.01

    # forward contango as a multiple of the current spot
    const forward_curve = [1.0, 1.025, 1.05, 1.075, 1.1, 1.125]
    const forward_contracts = 1:length(forward_curve)

    # Perishability
    const perishability = [1.0, 0.95, 0.9, 0.85, 0.8]
    const perisability_periods = 1:length(perishability)

    #==
        We need to discretise the price domain into "ribs". Our price varies
        between $3/widget and $9/widget, and we've chosen to have seven ribs.
        More ribs will give a more accurate solution, but result in more
        computation.
    ==#
    const ribs = collect(linspace(3, 9, 7))

    # A large number that we are confident is more profit that could be made
    # in the optimal solution
    const objecitve_bound = 10.0

end # @everywhere

#==
    Here is where we start the construction of our SDDP Model
==#
m = SDDPModel(
    # company is maximising profit
    sense             = :Max,
    # there are 12 months
    stages            = 12,
    # a large number to begin with
    objective_bound   = objecitve_bound,
    # a LP solver
    solver            = GurobiSolver(OutputFlag=0),
    #==
        We use Nested Average Value @ Risk for our risk measure.

            (1 - λ) * E[x] + λ * AV@R(1-β)[x]

        Increasing values of lambda are more risk averse.
        beta is the fraction of the tail we are worried about. Therefore,
        decreasing values of beta are more risk averse.
    ==#
    risk_measure      = NestedAVaR(lambda=0.5, beta=0.25),
    # price risk magic
    value_function    = InterpolatedValueFunction(
                            #==
                                Note: dynamics can't depend on other things
                                    and you must supply a function
                                    f(price::Float64,
                                        noise,
                                        stage::Int,
                                        markovstate::Int
                                    )
                                If you don't know what the markov state is, just
                                add the input, but don't reference it. It won't
                                matter.
                            ==#
                            dynamics       = price_dynamics,
                            # the initial price
                            initial_price  = initial_price,
                            # our rib discretisations
                            rib_locations  = ribs,
                            # and the noise for our price process
                            noise          = Noise(sampled_errors)
                        )
                                            #==
                                                so here 'sp' is a JuMP model
                                                for each subproblem, and 't' is
                                                an integer count from 1 to 12
                                                (the number of stages).
                                            ==#
                                            ) do sp, t

    # create state variables
    @states(sp, begin
        #==
            The number of contracts due 1, 2, 3, 4, 5, and 6 months from now.
            contacts0[1] is the number of contracts due for sale in this month).
        ==#
        contracts[month = forward_contracts]  >= 0, contracts0 == 0
        #==
            The total number of unsold widgets we have in storage at each
            perishablity level.
        ==#
        storage[p=perisability_periods] >= 0, storage0 == 0
    end)

    # auxillary variables
    @variables(sp, begin
        # Quantity of widgets to sell on spot
        sell_on_spot[p=perisability_periods]        >= 0
        # Quantity of product to deliver on contract
        deliver_to_contract[p=perisability_periods] >= 0

        # quantity of contracts to forward sell
        0 <= contract_sells[month = forward_contracts] <= 10
        # quantity to purchase on spot to satisfy contract
        buy_on_spot >= 0
        # production
        production >= 0
        # dummy variables
        widget_equivalent_sold
        total_contracts_traded
    end)

    # constraints
    @constraints(sp, begin
        # Each month we shift the contacts by 1 time period, and add any sells
        [i=forward_contracts[1:end-1]], contracts[i] ==
            contracts0[i+1] + contract_sells[i]
        # In the last period, it's just how many we sell this month
        contracts[forward_contracts[end]] ==
            contract_sells[forward_contracts[end]]

        # Each month unsold product perishes a bit more
        [i=perisability_periods[2:end]], storage[i] == storage0[i-1] -
            deliver_to_contract[i] - sell_on_spot[i]
        storage[1] == production - deliver_to_contract[1] - sell_on_spot[1] +
            buy_on_spot

        # We have to deliver the contract
        contracts[1] == sum(deliver_to_contract[p] for p in perisability_periods)

        # a dummy variable to make the objective simpler to recalculate
        widget_equivalent_sold ==
            sum(perishability[p] * sell_on_spot[p]
                for p in perisability_periods) -
            1.5 * buy_on_spot +
            sum(forward_curve[i] * contract_sells[i]
                for i in forward_contracts[end])

        total_contracts_traded == sum(contract_sells[i]
                                    for i in forward_contracts)

        # can't contract past end of year
        end_of_year[i=forward_contracts; i > 12-t], contract_sells[i] == 0
    end)

    # a constraint with varying RHS (but we leverage the JuMP tooling to evaluate that)
    @scenario(sp,
        alpha = linspace(0., 0.05, 10),
        production <= alpha
    )

    stageobjective!(sp, price -> price * widget_equivalent_sold -
        transaction_cost * total_contracts_traded
    )

end

SDDP.solve(m,
    # maximum number of cuts
    max_iterations = 100,
    # time limit for solver (seconds)
    time_limit     = 3600,
    # control the forward simulation part of the algorithm
    simulation = MonteCarloSimulation(
        # how often should we test convergence?
        frequency = 5,
        #==
            Simulation incrementation. Do min number of simulations, if there is
            evidence of convergence, do step more. Re-test. If there isn't any
            evidence of convergence keep going with the cutting iterations
            sice more simulations will only refine the estimate.

            If we reach max simulations and there is still evidence of
            convergence, then we can terminate with termination=true.
        ==#
        min         = 500,
        step        = 100,
        max         = 500,
        termination = false
    ),
    # if we have multiple processors loaded, use async solver
    solve_type      = nprocs()>2?Asyncronous():Serial(),
    reduce_memory_footprint = false,
    # write the log output to file
    log_file        = "contracting.log",
    # save all the discovered cuts to file as well
    # cut_output_file = "contracting.rib.cuts"
)

#==
    Simulate the policy 200 times
==#
results = simulate(m,
    200,               # number of monte carlo realisations
    [
        :contracts, :contracts0, :storage, :storage0, :sell_on_spot,
        :contract_sells, :deliver_to_contract, :buy_on_spot, :production, :price
    ]       # variables to return
)

#==
    Plot everything
==#
@visualise(results, i, t, begin
    results[i][:stageobjective][t], (title="Accumulated Profit",
                            ylabel="Accumulated Profit (\$)", cumulative=true)

    results[i][:price][t], (title="Price")

    results[i][:buy_on_spot][t], (title="Spot Buys")

    results[i][:production][t], (title="Production")

    results[i][:contracts0][t][1], (title="Contracts0 1")
    results[i][:contracts0][t][2], (title="Contracts0 2")
    results[i][:contracts0][t][3], (title="Contracts0 3")
    results[i][:contracts0][t][4], (title="Contracts0 4")
    results[i][:contracts0][t][5], (title="Contracts0 5")
    results[i][:contracts0][t][6], (title="Contracts0 6")

    results[i][:contracts][t][1], (title="Contracts 1")
    results[i][:contracts][t][2], (title="Contracts 2")
    results[i][:contracts][t][3], (title="Contracts 3")
    results[i][:contracts][t][4], (title="Contracts 4")
    results[i][:contracts][t][5], (title="Contracts 5")
    results[i][:contracts][t][6], (title="Contracts 6")

    results[i][:contract_sells][t][1], (title="Contract Sells 1")
    results[i][:contract_sells][t][2], (title="Contract Sells 2")
    results[i][:contract_sells][t][3], (title="Contract Sells 3")
    results[i][:contract_sells][t][4], (title="Contract Sells 4")
    results[i][:contract_sells][t][5], (title="Contract Sells 5")
    results[i][:contract_sells][t][6], (title="Contract Sells 6")

    results[i][:storage][t][1], (title="Storage 1")
    results[i][:storage][t][2], (title="Storage 2")
    results[i][:storage][t][3], (title="Storage 3")
    results[i][:storage][t][4], (title="Storage 4")
    results[i][:storage][t][5], (title="Storage 5")

    results[i][:storage0][t][1], (title="Storage0 1")
    results[i][:storage0][t][2], (title="Storage0 2")
    results[i][:storage0][t][3], (title="Storage0 3")
    results[i][:storage0][t][4], (title="Storage0 4")
    results[i][:storage0][t][5], (title="Storage0 5")

    results[i][:sell_on_spot][t][1], (title="Spot Sells 1")
    results[i][:sell_on_spot][t][2], (title="Spot Sells 2")
    results[i][:sell_on_spot][t][3], (title="Spot Sells 3")
    results[i][:sell_on_spot][t][4], (title="Spot Sells 4")
    results[i][:sell_on_spot][t][5], (title="Spot Sells 5")

    results[i][:deliver_to_contract][t][1], (title="Contract Delivery 1")
    results[i][:deliver_to_contract][t][2], (title="Contract Delivery 2")
    results[i][:deliver_to_contract][t][3], (title="Contract Delivery 3")
    results[i][:deliver_to_contract][t][4], (title="Contract Delivery 4")
    results[i][:deliver_to_contract][t][5], (title="Contract Delivery 5")
end)

#==
    Let's save the model to disk so we can come back to it later. Note this is
    pretty fragile and is not guaranteed to break between Julia serialiser
    versions.
==#
# SDDP.savemodel!(m, "contracting.sddpm")
#
# m2 = SDDP.loadsddpmodel("contracting.sddpm")
# results2 = simulate(m2,
#     200,               # number of monte carlo realisations
#     [
#         :contracts, :contracts0, :storage, :sell_on_spot,
#         :contract_sells, :buy_on_spot, :production
#     ]       # variables to return
# )
