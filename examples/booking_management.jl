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

function booking_management_model(num_days, num_rooms, num_requests, mip_solver)
    # maximum revenue that could be accrued.
    max_revenue = (num_rooms + num_requests) * num_days * num_rooms

    #=
        booking_requests is a vector of {0,1} arrays of size
        (num_days x num_rooms) if the room is requested.
    =#
    booking_requests = Array{Int, 2}[]
    for room in 1:num_rooms
        for day in 1:num_days
            # note: length_of_stay is 0 indexed to avoid unncecessary +/- 1
            # on the indexing
            for length_of_stay in 0:(num_days - day)
                req = zeros(Int, (num_rooms, num_days))
                req[room:room, day .+ (0:length_of_stay)] .= 1
                push!(booking_requests, req)
            end
        end
    end

    model = SDDP.LinearPolicyGraph(
            stages = num_requests, upper_bound = max_revenue,
            sense = :Max, optimizer = with_optimizer(GLPK.Optimizer),
            mip_solver = mip_solver
            ) do sp, stage

        @variable(sp,
            0 <= vacancy[room=1:num_rooms, day=1:num_days] <= 1,
            SDDP.State, Bin, initial_value = 1)
        @variables(sp, begin
            # Accept request for booking of room for length of time.
            0 <= accept_request <= 1, Bin
            # Accept a booking for an individual room on an individual day.
            0 <= room_request_accepted[1:num_rooms, 1:num_days] <= 1, Bin
            # Helper for JuMP.fix
            req[1:num_rooms, 1:num_days]
        end)
        for room in 1:num_rooms, day in 1:num_days
            @constraints(sp, begin
                # Update vacancy if we accept a room request
                vacancy[room, day].out == vacancy[room, day].in - room_request_accepted[room, day]
                # Can't accept a request of a filled room
                room_request_accepted[room, day] <= vacancy[room, day].in
                # Can't accept invididual room request if entire request is declined
                room_request_accepted[room, day] <= accept_request
                # Can't accept request if room not requested
                room_request_accepted[room, day] <= req[room, day]
                # Accept all individual rooms is entire request is accepted
                room_request_accepted[room, day] + (1-accept_request) >= req[room, day]
            end)
        end
        SDDP.parameterize(sp, booking_requests) do request
            JuMP.fix.(req, request)
        end
        @stageobjective(sp, sum(
            (room + stage - 1) * room_request_accepted[room, day]
            for room in 1:num_rooms for day in 1:num_days
            )
        )
    end
end

function booking_management(mip_solver)
    m_1_2_5 = booking_management_model(1, 2, 5, mip_solver)
    SDDP.train(m_1_2_5, iteration_limit = 10, print_level = 0)
    if mip_solver == SDDP.ContinuousRelaxation()
        @test SDDP.calculate_bound(m_1_2_5) >= 7.25
    else
        @test isapprox(SDDP.calculate_bound(m_1_2_5), 7.25, atol=0.02)
    end

    m_2_2_3 = booking_management_model(2, 2, 3, mip_solver)
    SDDP.train(m_2_2_3, iteration_limit = 40, print_level = 0)
    if mip_solver == SDDP.ContinuousRelaxation()
        @test SDDP.calculate_bound(m_1_2_5) > 6.13
    else
        @test isapprox(SDDP.calculate_bound(m_2_2_3), 6.13, atol=0.02)
    end
end

for mip_solver in [SDDP.SDDiP(with_optimizer(GLPK.Optimizer)), SDDP.ContinuousRelaxation()]
    booking_management(mip_solver)
end
