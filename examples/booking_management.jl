#  Copyright 2017-19, Oscar Dowson, Eyob Zewdie
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#=
    This example concerns the acceptance of booking requests for rooms in a
    hotel in the leadup to a large event.

    Each stage, we receive a booking request and can choose to accept or decline
    it. Once accepted, bookings cannot be terminated.
=#

using SDDP, GLPK, Test

function booking_management_model(NUMBER_OF_DAYS, NUMBER_OF_ROOMS, NUMBER_OF_REQUESTS)
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
                booking_request[room:room, day .+ (0:length_of_stay)] .= 1
                push!(BOOKING_REQUESTS, booking_request)
            end
        end
    end

    model = SDDP.LinearPolicyGraph(
            stages = NUMBER_OF_REQUESTS, upper_bound = MAX_REVENUE,
            sense = :Max, optimizer = with_optimizer(GLPK.Optimizer)
            ) do sp, stage
        @variable(sp,
            0 <= vacancy[room=1:NUMBER_OF_ROOMS, day=1:NUMBER_OF_DAYS] <= 1,
            SDDP.State, initial_value = 1)
        @variables(sp, begin
            # Accept request for booking of room for length of time.
            0 <= accept_request <= 1, Bin
            # Accept a booking for an individual room on an individual day.
            0 <= room_request_accepted[1:NUMBER_OF_ROOMS, 1:NUMBER_OF_DAYS] <= 1, Bin
            # Helper for JuMP.fix
            booking_request[1:NUMBER_OF_ROOMS, 1:NUMBER_OF_DAYS]
        end)
        for room in 1:NUMBER_OF_ROOMS, day in 1:NUMBER_OF_DAYS
            @constraints(sp, begin
                # Update vacancy if we accept a room request
                vacancy[room, day].out == vacancy[room, day].in - room_request_accepted[room, day]
                # Can't accept a request of a filled room
                room_request_accepted[room, day] <= vacancy[room, day].in
                # Can't accept invididual room request if entire request is declined
                room_request_accepted[room, day] <= accept_request
                # Can't accept request if room not requested
                room_request_accepted[room, day] <= booking_request[room, day]
                # Accept all individual rooms is entire request is accepted
                room_request_accepted[room, day] + (1-accept_request) >= booking_request[room, day]
            end)
        end
        SDDP.parameterize(sp, BOOKING_REQUESTS) do request
            JuMP.fix.(booking_request, request)
        end
        @stageobjective(sp, sum(
            (room + stage - 1) * room_request_accepted[room, day]
            for room in 1:NUMBER_OF_ROOMS for day in 1:NUMBER_OF_DAYS
            )
        )
    end
end

function booking_management()
    m_1_2_5 = booking_management_model(1, 2, 5)
    SDDP.train(m_1_2_5, iteration_limit = 10)
    @test isapprox(SDDP.calculate_bound(m_1_2_5), 7.25, atol=0.02)

    m_2_2_3 = booking_management_model(2, 2, 3)
    SDDP.train(m_2_2_3, iteration_limit = 40)
    @test isapprox(SDDP.calculate_bound(m_2_2_3), 6.13, atol=0.02)
end

booking_management()