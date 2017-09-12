#  Copyright 2017, Eyob Zewdie, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    This example concerns the acceptance of booking requests for rooms in a
    hotel in the leadup to a large event.

    Each stage, we receive a booking request and can choose to accept or decline
    it. Once accepted, bookings cannot be terminated.

    Properly formulated, this problem requires binary state variables. However,
    SDDP.jl cannot solve such problems. Therefore we solve the LP relaxation.

    Checkout github.com/lkapelevich/SDDiP.jl for a multi-stage stochastic
    programming solver with binary state variables.
=#
using SDDP, JuMP, Clp
using Base.Test

# Number of days for the event
const NUMBER_OF_DAYS      = 6
# Number of rooms in the hotel
const NUMBER_OF_ROOMS     = 4
# Number of booking requests we consider in the lead-up
const NUMBER_OF_REQUESTS  = 7

# revenue accrued if room `room` is booked on day `day`
const REVENUE = [
    room - 0.5 + (day-1)*0.2
    for room in 1:NUMBER_OF_ROOMS, day in 1:NUMBER_OF_DAYS
]

# maximum revenue that could be accrued.
const MAX_REVENUE = sum(REVENUE)

#=
    BOOKING_REQUESTS is a vector with one element for each stage (i.e.
    request). That element is a vector of {0,1} arrays of size
    (NUMBER_OF_DAYS x NUMBER_OF_ROOMS) if the room is requested.
=#
const BOOKING_REQUESTS = [
    Array{Int, 2}[]
    for request in 1:NUMBER_OF_REQUESTS
]
for request in 1:NUMBER_OF_REQUESTS
    for room in 1:NUMBER_OF_ROOMS
        for day in 1:NUMBER_OF_DAYS
            # note: length_of_stay is 0 indexed to avoid unncecessary +/- 1
            # on the indexing
            for length_of_stay in 0:(NUMBER_OF_DAYS - day)
                # ==========================================================
                #   put some stay length logic here
                if length_of_stay + 1  < room
                    # min of 1 nights in room 1, 2 in room 2, etc
                    continue
                end
                # ==========================================================
                booking_request = zeros(Int, (NUMBER_OF_ROOMS, NUMBER_OF_DAYS))
                booking_request[room:room, day + (0:length_of_stay)] = 1
                push!(BOOKING_REQUESTS[request], booking_request)
            end
        end
    end
end

m = SDDPModel(
             stages = NUMBER_OF_REQUESTS,
    objective_bound = MAX_REVENUE,
              sense = :Max,
             solver = ClpSolver()
                        ) do sp, stage

    # Seat remaining?
    @state(sp, 0 <= vacancy_after_decision[room=1:NUMBER_OF_ROOMS, day=1:NUMBER_OF_DAYS] <= 1, vacancy_before_decision==1)

    @variables(sp, begin
        # accept request for booking of room for length of time
        0 <= accept_request <= 1
        # accept a booking for an individual room on an individual day
        0 <= room_request_accepted[1:NUMBER_OF_ROOMS, 1:NUMBER_OF_DAYS] <= 1
    end)

    for room in 1:NUMBER_OF_ROOMS, day in 1:NUMBER_OF_DAYS
        @constraints(sp, begin
            # Update vacancy if we accept a room request
            vacancy_after_decision[room, day] == vacancy_before_decision[room, day] - room_request_accepted[room, day]
            # can't accept a request of a filled room
            room_request_accepted[room, day] <= vacancy_before_decision[room, day]
            # can't accept invididual room request if entire request is declined
            room_request_accepted[room, day] <= accept_request

        end)
        @noises(sp, booking_request = BOOKING_REQUESTS[stage], begin
            # can't accept request if room not requested
            room_request_accepted[room, day] <= booking_request[room, day]
            # accept all individual rooms is entire request is accepted
            room_request_accepted[room, day] + (1-accept_request) >= booking_request[room, day]
        end)
    end

    stageobjective!(sp, sum(
        REVENUE[room, day]*room_request_accepted[room, day]
        for room in 1:NUMBER_OF_ROOMS for day in 1:NUMBER_OF_DAYS
        )
    )

end

solution = solve(m,
     max_iterations = 75,
    #  print_level = 0,
     # simulate the lower bound of the policy
     simulation = MonteCarloSimulation(
                    frequency = 50,
                          min = 100,
                         step = 100,
                          max = 500
      )
)

srand(111)

# I get different answers running this same file with Clp
@test isapprox(getbound(m), 37.0, atol=5.0)

states = Union{Float64, Vector{Float64}}[0.0 for i in 1:NUMBER_OF_ROOMS, j in 1:NUMBER_OF_DAYS]
# compare (room, day) = (1,1)
states[1,1] = collect(linspace(0, 1, 51))
# against (room, day) = (2,1)
states[2,1] = collect(linspace(0, 1, 51))

# ==============================================================================
#   Uncomment to plot value function
# SDDP.plotvaluefunction(m, 1, 1, states..., label1 = "Room (1,1)", label2 = "Room (2,1)")
