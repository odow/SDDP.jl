#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

# Add procs for parallel speeeeeeed
addprocs(3)

# only need this on master
using Base.Test

@everywhere begin
    # need these everywhere
    using SDDP, JuMP, Clp

    srand(myid() * 11111)

    # Demand for newspapers
    # There are two equally probable noises in each stage
    #   Demand[stage, noise]
    Demand = [
        10.0 15.0;
        12.0 20.0;
         8.0 20.0
    ]

    # Markov state purchase prices
    PurchasePrice = [5.0, 8.0]

    RetailPrice = 7.0

    # Transition matrix
    Transition = Array{Float64, 2}[
        [1.0]',
        [0.6 0.4],
        [0.3 0.7; 0.3 0.7]
      ]

      # overload to check parallel
      function SDDP.storekey!(::Type{Val{:pid}}, store, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int)
          push!(store, myid())
      end
end

function createmodel(risk_measure)
    # Initialise SDDP Model
    m = SDDPModel(
            sense             = :Max,
            stages            = 3,
            objective_bound   = 1000,
            markov_transition = Transition,
            risk_measure      = risk_measure,
            solver          = ClpSolver()
                                                    ) do sp, stage, markov_state

        # ====================
        #   State variable
        @state(sp, 0 <= stock <= 100, stock0==5)

        # ====================
        #   Other variables
        @variables(sp, begin
            buy  >= 0  # Quantity to buy
            sell >= 0  # Quantity to sell
        end)

        # ====================
        #   Noises
        @rhsnoises(sp, D=Demand[stage,:], begin
            sell <= D
            sell >= 0.5D
        end)

        # ====================
        #   Objective
        @stageobjective(sp, sell * RetailPrice - buy * PurchasePrice[markov_state])

        # ====================
        #   Dynamics constraint
        @constraint(sp, stock == stock0 + buy - sell)

    end
end

m = createmodel(NestedAVaR(beta   = 0.6, lambda = 0.5))


@test_throws Exception SDDP.solve(m, max_iterations=30, solve_type=Asyncronous(slaves=[2]))

# slave processes
slaves = workers()
#=
 Add 1 slave after every 2 iterations until all are used. i.e.
    Iteration | 1 | 2 | 3 | 4 | 5 | 6 | 7 | ...
    # Slaves  | 1 | 1 | 2 | 2 | 3 | 3 | 4 | ...
=#
slave_steps = 2.0
solvestatus = SDDP.solve(m,
    max_iterations = 30,
    solve_type     = Asyncronous(slaves=slaves, step=slave_steps),
   cut_output_file = "async.cuts",
    simulation     = MonteCarloSimulation(
                        frequency = 10,
                        min       = 5,
                        max       = 50,
                        step      = 5
                             )
)

@test isapprox(getbound(m), 93.267, atol=1e-3)
@test solvestatus == :max_iterations

sim = simulate(m, 100, [:stock, :pid])
@test length(sim) == 100

# check pid's are not all 1 so that some were simulated on other cores
pids = [s[:pid][1] for s in sim]
@test !all(pids .== 1)

m4 = createmodel(NestedAVaR(beta   = 0.6, lambda = 0.5))
loadcuts!(m4, "async.cuts")
rm("async.cuts")

SDDP.solve(m4, max_iterations=1, solve_type=Asyncronous())
@test isapprox(getbound(m), getbound(m4), atol=1e-3)

m2 = createmodel(Expectation())

solvestatus = SDDP.solve(m2,
    max_iterations = 50, print_level=0,
    solve_type     = Asyncronous(),
    cut_selection_frequency = 5,
    simulation     = MonteCarloSimulation(
                        frequency = 10,
                        min       = 5,
                        max       = 50,
                        step      = 5,
                        termination = true
                             )
)

@test solvestatus == :converged


m3 = createmodel(Expectation())

solvestatus = SDDP.solve(m3,
    time_limit = 0.1, print_level=0,
    solve_type     = Asyncronous(slaves=vcat(workers(), myid()))
)

@test solvestatus == :time_limit

# on Julia v0.5 waitfor defaults to 0.0 ...
rmprocs(workers(), waitfor=60.0)
