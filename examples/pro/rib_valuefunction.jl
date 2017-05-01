#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp

# For repeatability
srand(11111)

sampled_errors = [-0.1290, -0.1010, -0.0814, -0.0661, -0.0530, -0.0412, -0.0303, -0.0199, -0.00987, 0.0, 0.00987, 0.0199, 0.0303, 0.0412, 0.0530, 0.0661, 0.0814, 0.1010, 0.1290]

σ² = linspace(1, 0, 28) # Decreasing variance in changes in price over time
transaction_cost = 0.01

markov_states = [0.9, 1.1]

markov_transition = Array{Float64, 2}[
    [0.5 0.5]
]
for t in 2:28
    push!(markov_transition, [0.75 0.25; 0.3 0.7])
end

box(x, a, b) = min(b, max(a, x))
price_dynamics(p, w, t) = box(exp(log(p) + σ²[t]*w), 3, 9)

m = SDDPModel(
    stages       = 28,
    markovstates = 2,
    transition = markov_transition,
    scenarios = ones(19) / 19,
    sense        = :Max
                        ) do sp, t, i

    # create state variables
    @states(sp, begin
        0 <= contracts  <= 1.5, contracts0 == 0
        0 <= production <= 1.5, production0 == 0
    end)

    # auxillary variables
    @variables(sp, begin
        0 <= buy <= 1.2
        0 <= sell <= 1.2
        output >= 0
    end)

    # constraints
    @constraints(sp, begin
        contracts  == contracts0 + buy - sell
        production == production0 + output
    end)

    # a constraint with varying RHS (but we leverage the JuMP tooling to evaluate that)
    @scenario(sp,
        alpha = sampled_errors,
        output <= alpha * markov_states[i]
    )

    if t < 28
        stageobjective!(sp,
            linspace(3, 9, 5),                  # locations
            (p, w)->price_dynamics(p, w, t),    # dynamics for updating price state
            linspace(0, 0.05, 5),   # noises
            fill(0.2, 5),           # probability support
            price -> (buy * price - transaction_cost * (buy + sell)) # returns AffExpr for stage objective
        )
    else
        stageobjective!(sp,
            linspace(3, 9, 5),                  # locations
            (p, w)->price_dynamics(p, w, t),    # dynamics for updating price state
            linspace(0, 0.05, 5), # noises
            # fill(0.2, 5),           # probability support
            price -> (production - contracts) * price # objective
        )
    end
end
