#  Copyright 2017,  Oscar Dowson, Eyob Zewdie
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
This problem is a version of the Ambulance dispatch problem. A hospital is
located at 0 on the number line that stretches from 0 to 100. Ambulance bases
are located at points 20, 40, 60, 80, and 100. When not responding to a call,
Ambulances must be located at a base, or the hospital. In this example there are
three ambulances.

Example location:

    H       B       B       B       B       B
    0 ---- 20 ---- 40 ---- 60 ---- 80 ---- 100

Each stage, a call comes in from somewhere on the number line. The agent must
decide which ambulance to dispatch. They pay the cost of twice the driving
distance. If an ambulance is not dispatched in a stage, the ambulance can be
relocated to a different base in preparation for future calls. This incurrs a
cost of the driving distance.

Properly formulated, this problem requires binary state variables. However,
SDDP.jl cannot solve such problems. Therefore we solve the LP relaxation.

Checkout github.com/lkapelevich/SDDiP.jl for a multi-stage stochastic
programming solver with binary state variables.
=#
using SDDP, JuMP, Clp, Base.Test

function vehiclelocationmodel(nvehicles, baselocations, requestlocations)
    Locations = baselocations # Points on the number line where emergency bases are located
    Vehicles = 1:nvehicles          # ambulances
    Bases = 1:length(baselocations) # base 1 = hostpital
    Requests = collect(0:2:100)     # Points on the number line where calls come from

    shiftcost(src, dest) = abs(Locations[src] - Locations[dest])
    dispatchcost(request, base) = 2 * (abs(request - Locations[1]) + abs(request-Locations[base]))

    #Initial State of emergency vehicles at bases
    Q0 = zeros(Int64, (length(Bases), length(Vehicles)))
    Q0[1,:] = 1 # all ambulances start at hospital

    m = SDDPModel(
                 stages = 10,
        objective_bound = 0.0,
                  sense = :Min,
                 solver = ClpSolver()
                            ) do sp, t

        # Vehicles at which bases?
        @state(sp, 0 <= q[b=Bases, v=Vehicles] <= 1, q0 == Q0[b, v])

        @variables(sp, begin
            # which vehicle is dipatched?
            0 <= dispatch[src=Bases, v=Vehicles] <= 1
            # shifting vehicles between bases
            0 <= shift[src=Bases, v=Vehicles, dest=Bases] <= 1
        end)

        @expression(sp, basebalance[b in Bases, v in Vehicles],
            # flow of vehicles in and out of bases:
            #initial - dispatch      - shifted from      + shifted to
            q0[b, v] - dispatch[b,v] - sum(shift[b,v,:]) + sum(shift[:,v,b])
        )
        @constraints(sp, begin
            # only one vehicle dispatched to call
            sum(dispatch) == 1
            # can only dispatch vehicle from base if vehicle is at that base
            [b in Bases, v in Vehicles], dispatch[b,v] <= q0[b,v]
            # can only shift vehicle if vehicle is at that src base.
            [b in Bases, v in Vehicles], sum(shift[b,v,:]) <= q0[b, v]
            # can only shift vehicle if vehicle is not being dispatched
            [b in Bases, v in Vehicles], sum(shift[b,v,:]) + dispatch[b,v] <= 1
            # can't shift to same base
            [b in Bases, v in Vehicles], shift[b,v,b] == 0
            # Update states for non-home/non-hospital bases
            [b in Bases[2:end], v in Vehicles], q[b, v] == basebalance[b,v]
            # Update states for home/hospital bases
            [v in Vehicles], q[1, v] == basebalance[1,v] + sum(dispatch[:,v])
        end)

        @stageobjective(sp, request=Requests, sum(
                #distance to travel from base to emergency and back to home base
                dispatch[b,v] * dispatchcost(request, b) +
                #distance travelled by vehilces relocating bases
                sum(shiftcost(b, dest) * shift[b, v, dest] for dest in Bases)
            for b in Bases, v in Vehicles)
        )
    end
end

ambulancemodel = vehiclelocationmodel(3, [0, 20, 40, 60, 80, 100], 0:2:100)
@test solve(ambulancemodel, max_iterations=50) == :max_iterations
@test isapprox(getbound(ambulancemodel), 1700.0, atol=5)
# N = 100
# results = simulate(m, N, [:dispatch, :q, :shift])
#
# plt = SDDP.newplot()
# for b in Bases
#     SDDP.addplot!(plt, 1:N, 1:T, (i,t)->sum(results[i][:q][t][b,:]), title="Base $(b)")
# end
# SDDP.addplot!(plt, 1:N, 1:T, (i,t)->sum(results[i][:shift][t]), title="Shifts")
# SDDP.addplot!(plt, 1:N, 1:T, (i,t)->Requests[results[i][:noise][t]], title="Request")
# SDDP.show(plt)
