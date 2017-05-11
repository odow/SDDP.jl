#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==

Example: The Widget Producer

A company produces divisible widgets (i.e. a continuous product, rather than
discrete units) for sale on a spot market each month. The quantity of widgets
they produce is uncertain and comes from some finite distribution.

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
using SDDP, JuMP, Clp, Base.Test

# For repeatability
srand(11111)

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

#==
    There is a transaction cost to trading in the forward contracts
==#
const transaction_cost = 0.01

#==
    We need to discretise the price domain into "ribs". Our price varies between
    $3/widget and $/widget, and we've chosen to have five ribs. More ribs will
    give a more accurate solution, but result in more computation.
==#
ribs = collect(linspace(3, 9, 5))



m = SDDPModel(
    # company is maximising profit
    sense             = :Max,
    # there are 12 months
    stages            = 12,
    # a large number to begin with
    objective_bound   = 1e9,
    # a LP solver
    solver            = ClpSolver(),
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
            contacts0 is the number of contracts due at the end of the last
            stage (so that contracts[1] is the number of widgets contracted for
            sale in this month).
        ==#
        0 <= contracts[month = 1:6]  <= 1.5, contracts0 == 0
        #==
            The total number of unsold widgets we have in storage
        ==#
        0 <= production <= 1.5, production0 == 0
    end)

    # auxillary variables
    @variables(sp, begin
        # quantity of widgets to sell on spot
        sell_on_spot >= 0
        # 
        output >= 0
    end)

    # constraints
    @constraints(sp, begin
        contracts_sold  == contracts_sold0 + sell
        production == production0 + output
        contracts_sold0 + sell <= 1.2
    end)

    # a constraint with varying RHS (but we leverage the JuMP tooling to evaluate that)
    @scenario(sp,
        alpha = linspace(0., 0.05, 5),
        output <= alpha
    )

    if t < 28
        stageobjective!(sp,
            price -> (sell * (price - transaction_cost)) # returns AffExpr for stage objective
        )
    else
        stageobjective!(sp,
            price -> (production - contracts_sold0) * price # objective
        )
    end
end

SDDP.solve(m,
    max_iterations = 5,
    simulation = MonteCarloSimulation(
        frequency = 5,
        min       = 500,
        max       = 500
    )
)

@test getbound(m) >= 4.1
@test getbound(m) <= 4.3
