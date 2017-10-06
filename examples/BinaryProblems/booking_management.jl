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
using SDDP, JuMP, Clp, Base.Test

function bookingmanagementmodel(NUMBER_OF_DAYS, NUMBER_OF_ROOMS, NUMBER_OF_REQUESTS)

    # maximum revenue that could be accrued.
    MAX_REVENUE = (NUMBER_OF_ROOMS + NUMBER_OF_REQUESTS) * NUMBER_OF_DAYS * NUMBER_OF_ROOMS

    #=
        BOOKING_REQUESTS is a vector of {0,1} arrays of size
        (NUMBER_OF_DAYS x NUMBER_OF_ROOMS) if the room is requested.
    =#
    BOOKING_REQUESTS = Array{Int, 2}[]
    for room in 1:NUMBER_OF_ROOMS
        for day in 1:NUMBER_OF_DAYS
            # note: length_of_stay is 0 indexed to avoid unncecessary +/- 1
            # on the indexing
            for length_of_stay in 0:(NUMBER_OF_DAYS - day)
                booking_request = zeros(Int, (NUMBER_OF_ROOMS, NUMBER_OF_DAYS))
                booking_request[room:room, day + (0:length_of_stay)] = 1
                push!(BOOKING_REQUESTS, booking_request)
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
            @rhsnoises(sp, booking_request = BOOKING_REQUESTS, begin
                # can't accept request if room not requested
                room_request_accepted[room, day] <= booking_request[room, day]
                # accept all individual rooms is entire request is accepted
                room_request_accepted[room, day] + (1-accept_request) >= booking_request[room, day]
            end)
        end
        @stageobjective(sp, sum(
            (room + stage - 1) * room_request_accepted[room, day]
            for room in 1:NUMBER_OF_ROOMS for day in 1:NUMBER_OF_DAYS
            )
        )
    end
end

srand(1234)
m_1_2_5 = bookingmanagementmodel(1, 2, 5)
@test solve(m_1_2_5, max_iterations = 10, print_level = 0) == :max_iterations
@test isapprox(getbound(m_1_2_5), 7.25, atol=0.001)

srand(1234)
m_2_2_3 = bookingmanagementmodel(2, 2, 3)
@test solve(m_2_2_3, max_iterations = 40, print_level = 2) == :max_iterations
@test isapprox(getbound(m_2_2_3), 6.13, atol=0.001)
