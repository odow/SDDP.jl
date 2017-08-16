#Managing booking requests for rooms in a hotel
using SDDP, JuMP, Clp

srand(111111)

const NUMBER_OF_DAYS      = 6 # Number of days
const NUMBER_OF_ROOMS     = 4 # Number of rooms (room 1 cost 1 dollar, room 2 cost 2 dollar etc.)
const NUMBER_OF_REQUESTS  = 7 # Number of requests (i.e. stages)

# revenue
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
     max_iterations = 1000,
        # simulate the lower bound of the policy
         simulation = MonteCarloSimulation(
                        frequency = 50,
                              min = 100,
                             step = 100,
                              max = 500
              ),
)
