#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==

Example: The Perishable Widget Producer

A company produces perishable, divisible widgets (i.e. a continuous product,
rather than discrete units) for sale on a spot market each month. The quantity
of widgets they produce is uncertain and comes from some finite distribution.

The company can store the widgets, but, over time, the value of the widgets
decreases. Eventually the widgets perish and are worthless.

The spot price is determined by an auction system, and so varies from month to
month, but demonstrates serial correlation. In each auction, there is sufficient
demand that the widget producer finds a buyer for all their widgets, regardless
of the quantity they supply. Furthermore, the spot price is independent of the
widget producer (they are a small player in the market).

The spot price is highly volatile, and is the result of a process that is out of
the control of the company. To counteract their price risk, the company engages
in a forward contracting programme.

The forward contracting programme is a deal for physical widgets at a future
date in time. The company can engage in contracts for sales up to six months in
the future.

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

To run this example, there are two options.

1. From the command line:

    cmd> julia -p N perishablewidgets.jl

    This will start Julia with N processors. It N > 2, SDDP.jl will choose to
    solve using the Asynchronous solver. Upon completion, Julia will exit so make
    sure you write to file (or save) any results.

2. Interactively:
    cmd> julia -p N
    julia> include("perishablewidgets.jl")

    Upon completion, the model and results will be available for inspection.
==#

# These are a few of my favourite things...
using SDDP, JuMP, Gurobi

#==
    For repeatability we should set the random number seed. However, if we
    choose the same seed on all processors, then we end up doing identical
    things! myid() returns an integer of the pocessor id, so this way we get
    a unique random seed for all processors.
==#
@everywhere srand(myid() * 10)

#==
    Here is where we start the construction of our SDDP Model
==#
m = SDDPModel(
    # company is maximising profit
    sense             = :Max,
    # there are 12 months
    stages            = 12,
    # A large number that we are confident is more profit that could be made
    # in the optimal solution
    objective_bound   = 10.0,
    # a LP solver
    solver            = GurobiSolver(OutputFlag=0),
    # solver            = ClpSolver(),
    #==
        We use Nested Average Value @ Risk for our risk measure:

            (1 - λ) * E[X] + λ * AV@R(β)[X].

        Increasing values of lambda are more risk averse. beta is the fraction
        of the tail we are worried about. Therefore, decreasing values of beta
        are more risk averse.
    ==#
    risk_measure      = EAVaR(lambda=0.5, beta=0.25),
    # price risk magic
    value_function    = StaticPriceInterpolation(
                            #==
                                The spot price is a mean reverting, log-normal,
                                AR(1) model:

                                log(y(t)) = (1 - α) log(y(t-1)) + α * log(μ) + e,

                                where α is the rate at which the process reverts
                                to the mean, and μ is the mean that the process
                                reverts to.

                                Note: dynamics can't depend on other things and
                                you must supply a function
                                    f(price::Float64,
                                        noise,
                                        stage::Int,
                                        markovstate::Int
                                    )
                                If you don't know what the markov state is, just
                                add the input, but don't reference it. It won't
                                matter.
                            ==#
                            dynamics = (price, noise) -> begin
                                reversion_rate = 0.05
                                process_mean = 6.0
                                lower = 3.0
                                upper = 9.0
                                new_price = exp(
                                            (1 - reversion_rate) * log(price) +
                                            reversion_rate * log(process_mean) +
                                            noise
                                        )
                                return min(upper, max(lower, new_price))
                            end,
                            # The spot price in period 0
                            initial_price  = 4.50,
                            #==
                                We need to discretise the price domain into
                                "ribs". Our price varies between $3/widget and
                                $9/widget, and we've chosen to have seven ribs.
                                More ribs will give a more accurate solution,
                                but result in more computation.
                            ==#
                            rib_locations  =  collect(linspace(3, 9, 7)),
                            # and the noise for our price process
                            noise          = DiscreteDistribution([-0.1290, -0.1010, -0.0814,
                                                    -0.0661, -0.0530, -0.0412,
                                                    -0.0303, -0.0199, -0.00987,
                                                    0.0, 0.00987, 0.0199,
                                                    0.0303, 0.0412, 0.0530,
                                                    0.0661, 0.0814, 0.1010,
                                                0.1290]),
                            # and include some cut selection
                            cut_oracle     = LevelOneCutOracle()
                        )
                                            #==
                                                so here 'sp' is a JuMP model
                                                for each subproblem, and 't' is
                                                an integer count from 1 to 12
                                                (the number of stages).
                                            ==#
                                            ) do sp, t
    #==
        We begin by creating initialising some data.
        We're wrapping this inside the SDDPModel so that when we save it the
        objects get saved as well. Otherwise they just hold references to global
        objects and we'll need to initialise them each time as well.
    ==#

    #==
        There is a transaction cost to trading in the forward contracts

        A transaction cost also stops the model arbitraging numerical precision
        errors between the expected value of the S.A.A and the current price.
    ==#
    transaction_cost = 0.01

    # forward contango as a multiple of the current spot
    forward_curve = [1.0, 1.025, 1.05, 1.075, 1.1, 1.125]
    forward_contracts = 1:length(forward_curve)

    # Perishability
    perishability = [1.0, 0.95, 0.9, 0.85, 0.8]
    perishability_periods = 1:length(perishability)
    # initial stockpiles of product
    initial_storage = fill(0, length(perishability))

    #==
        Stochastic Production outcomes. Quantity of widgest produced in a stage.
    ==#
    production_realisations = linspace(0.1, 0.2, 5)

    #==
        create state variables. These need have the syntax

            outgoingstate, incomingstate == initial_value

        outgoingstate can be any valid JuMP @variable syntax including things
        such as
            0 <= state[i=1:3, j=1:2; i<j] <= 1

        incomingstate must be a single symbol, but the initial_value can use any
        iterators created by the outgoingstate variable (see storage0 definition)
    ==#
    @states(sp, begin
        #==
            The number of contracts due 1, 2, 3, 4, 5, and 6 months from now.
            contacts0[1] is the number of contracts due for sale in this month).
        ==#
        contracts[month = forward_contracts] >= 0, contracts0 == 0
        #==
            The total number of unsold widgets we have in storage at each
            perishablity level.
        ==#
        storage[p=perishability_periods] >= 0, storage0 == initial_storage[p]
    end)

    # any old JuMP variables
    @variables(sp, begin
        # Quantity of widgets to sell on spot
        sell_on_spot[p=perishability_periods]        >= 0
        # Quantity of product to deliver on contract
        deliver_to_contract[p=perishability_periods] >= 0

        #==
            Quantity of contracts to forward sell
            Put an upper bound on this to avoid the model arbitraging to Inf
        ==#
        0 <= contract_sells[month = forward_contracts] <= 10
        # quantity to purchase on spot to satisfy contract
        buy_on_spot >= 0
        # production
        production >= 0
        # dummy variables
        widget_equivalent_sold
        total_contracts_traded
    end)

    # any old JuMP constraints
    @constraints(sp, begin
        #==
            Each month we shift the contacts by 1 time period, and add any sells
            Note: this doesn't include contracts0[1] since those are the
            contracts to be delivered this month
        ==#
        [i=forward_contracts[1:end-1]], contracts[i] ==
            contracts0[i+1] + contract_sells[i]

        #==
            In the last period, it's just how many we sell this month
        ==#
        contracts[forward_contracts[end]] ==
            contract_sells[forward_contracts[end]]

        #==
            Each month unsold product perishes a bit more
        ==#
        [i=perishability_periods[2:end]], storage[i] == storage0[i-1] -
            deliver_to_contract[i] - sell_on_spot[i]
        #==
            But product made this period doesn't perish. This includes the
            production, less that which we delivered as fresh, less sold as
            fresh. We also assume that we buy fresh product on spot.
        ==#
        storage[1] == production - deliver_to_contract[1] - sell_on_spot[1] +
            buy_on_spot

        #==
            We have to deliver the contracts we agreed at the end of the last
            period. We can deliver from any perishability state, but we'll earn
            less (in the objective) for supplying more perished products.
        ==#
        contracts0[1] == sum(deliver_to_contract[p] for p in perishability_periods)

        #==
            A dummy variable to make the objective simpler to recalculate.
        ==#
        widget_equivalent_sold ==
            # we earn the perished product sold on spot
            sum(perishability[p] * sell_on_spot[p] -
                # plus the penalty on perished items sold on contract
                (1 - perishability[p]) * deliver_to_contract[p]
                for p in perishability_periods) -
            #==
                buying on spot incurrs some additional cost over and above
                the actual spot price (i.e. due to faster shipment times)
            ==#
            1.5 * buy_on_spot +
            # and we earn the forward price of each contract
            sum(forward_curve[i] * contract_sells[i]
                for i in forward_contracts)

        #==
            We use this to sum up transaction costs. This could go in the
            objective but then it'd make our stageobjective more complicated
            (although still linear).
        ==#
        total_contracts_traded == sum(contract_sells[i]
                                    for i in forward_contracts)

        # can't contract past end of year
        end_of_year[i=forward_contracts; i > 12-t], contract_sells[i] == 0
    end)

    #==
        Production of widgets is limited by some stochastic process
    ==#
    @rhsnoise(sp, alpha = production_realisations, production <= alpha)

    #==
        Set the stage objective. This should take function that maps the price
        into an AffExpr using the variables defined in the subproblem above.


    ==#
    @stageobjective(sp, price -> price * widget_equivalent_sold -
        transaction_cost * total_contracts_traded
    )

end

@time status = SDDP.solve(m,
    # maximum number of cuts to be added to any individual subproblem
    iteration_limit = 300,

    # time limit for solver (seconds)
    time_limit     = 20,

    #==
        Instead of the typical forwarad/backward pass split of the algorithm,
        I split it into cutting and simulation phases. In the cutting phase, we
        do forward/backward pass pairs (one iteration). Each iteration generates
        one cut at each subproblem.

        When we wish to test for convergence we switch to the simulation phase.
        This samples noises that have not, and are not, used in the cutting
        phase. The simulation keyword controls the simulation phase of the
        algorithm
    ==#
    simulation = MonteCarloSimulation(
        # how often should we test convergence?
        frequency = 100,
        #==
            Simulation incrementation. Do min number of simulations, if there is
            evidence of convergence, do step more. Re-test. If there isn't any
            evidence of convergence keep going with the cutting iterations
            sice more simulations will only refine the estimate.

            If we reach max simulations and there is still evidence of
            convergence, then we can terminate with terminate=true.
        ==#
        min       = 500,
        step      = 100,
        max       = 500,
        terminate = false
    ),

    #  SDDP.jl will automatically chose the best solver but we can orverride
    # solve_type      = Serial(),
    # solve_type      = Asynchronous(),

    # write the log output to file
    log_file        = "perishablewidgets.log",

    #==
        Save all the discovered cuts to file as well. This can be useful if you
        wish to analyse the cuts that are generated, or rebuild the model using
        them later. Another option is to use the SDDP.savemodel! function
        (discussed below).
    ==#
    # cut_output_file = "contracting.rib.cuts"

    #==
        Frequency with which the problems are rebuilt. No good rule. Aggressive
        cut selection reduces the size of the subproblems (good) but has
        additional overhead. However, a high cut_selection_frequency allows the
        subproblems to grow large which causes a big slow-down.
    ==#
    cut_selection_frequency = 20
)

#==
    Simulate the policy 200 times.

    note: :price is a special variable that is only valid for the
        InterpolatedValueFunction.

    results is a vector containing a dictionary for each simulation. In addition
    to the variables specified in the function call, other special keys are
        :stageobjective - costs incurred during the stage (not future)
        :obj            - :stageobjective + future value
        :markov         - index of markov state visited
        :noise       - index of noise visited
        :objective      - Total objective of simulation

    All values can be accessed as follows
        results[simulation index][key][stage]

    with the exception of :objective which is just
        results[simulation index][:objective]
==#
results = simulate(m,
    200,    # number of monte carlo realisations
    [       # variables to return
        :contracts, :contracts0, :storage, :storage0, :sell_on_spot,
        :contract_sells, :deliver_to_contract, :buy_on_spot, :production, :price
    ]
)

#==
    Plot everything using interactive javascript. This will launch an HTML page
    to explore.

    After the plot definition, keyword arguments can be used (in parenthesises):
        title       - set the title of the plot
        ylabel      - set the yaxis label
        xlabel      - set the xaxis label
        interpolate - interpolate lines between stages. Defaults to "linear"
        see https://github.com/d3/d3-3.x-api-reference/blob/master/SVG-Shapes.md
            #line_interpolate for all options
==#
plt = SDDP.newplot()
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:stageobjective][t],
  title="Accumulated Profit", ylabel="Accumulated Profit (\$)", cumulative=true)
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:price][t], title="Price")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:buy_on_spot][t], title="Spot Buys")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:production][t], title="Production")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts0][t][1], title="Contracts0 1")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts0][t][2], title="Contracts0 2")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts0][t][3], title="Contracts0 3")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts0][t][4], title="Contracts0 4")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts0][t][5], title="Contracts0 5")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts0][t][6], title="Contracts0 6")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts][t][1], title="Contracts 1")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts][t][2], title="Contracts 2")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts][t][3], title="Contracts 3")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts][t][4], title="Contracts 4")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts][t][5], title="Contracts 5")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contracts][t][6], title="Contracts 6")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contract_sells][t][1], title="Contract Sells 1")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contract_sells][t][2], title="Contract Sells 2")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contract_sells][t][3], title="Contract Sells 3")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contract_sells][t][4], title="Contract Sells 4")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contract_sells][t][5], title="Contract Sells 5")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:contract_sells][t][6], title="Contract Sells 6")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage][t][1], title="Storage 1")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage][t][2], title="Storage 2")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage][t][3], title="Storage 3")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage][t][4], title="Storage 4")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage][t][5], title="Storage 5")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage0][t][1], title="Storage0 1")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage0][t][2], title="Storage0 2")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage0][t][3], title="Storage0 3")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage0][t][4], title="Storage0 4")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:storage0][t][5], title="Storage0 5")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:sell_on_spot][t][1], title="Spot Sells 1")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:sell_on_spot][t][2], title="Spot Sells 2")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:sell_on_spot][t][3], title="Spot Sells 3")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:sell_on_spot][t][4], title="Spot Sells 4")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:sell_on_spot][t][5], title="Spot Sells 5")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:deliver_to_contract][t][1], title="Contract Delivery 1")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:deliver_to_contract][t][2], title="Contract Delivery 2")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:deliver_to_contract][t][3], title="Contract Delivery 3")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:deliver_to_contract][t][4], title="Contract Delivery 4")
SDDP.addplot!(plt, 1:200, 1:12, (i,t)->results[i][:deliver_to_contract][t][5], title="Contract Delivery 5")

#==
    Let's save the model to disk so we can come back to it later. Note this is
    pretty fragile and is not guaranteed to break between Julia serialiser
    versions.
==#
SDDP.savemodel!("contracting.sddpm", m)

# load a saved model
m2 = SDDP.loadmodel("contracting.sddpm")

# simulate with the saved model
results2 = simulate(m2,
    200,               # number of monte carlo realisations
    [   :contracts, :contracts0, :storage, :sell_on_spot,
        :contract_sells, :buy_on_spot, :production         ]
)

rm("contracting.sddpm")
rm("perishablewidgets.log")
