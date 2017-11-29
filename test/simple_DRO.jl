# #  Copyright 2017, Lea Kapelevich
#
# using SDDPPro, JuMP, Clp, Base.Test, SDDP
#
# # Look for the probability distribution when future costs are [1, 2, 3, 4, 5]
#
# const N = 5
#
# function buildandsolvemodel(measure::SDDP.AbstractRiskMeasure)
#
#     m = SDDPModel(
#     stages          = 2,
#     objective_bound = 1000,
#     sense           = :Max,
#     risk_measure    = measure,
#     solver          = ClpSolver()
#     ) do sp, stage
#
#         # ====================
#         #   State variable
#         @state(sp, x >= 0, x0 == 0)
#
#         # ====================
#         #   Uncertain variable
#         @variable(sp, ω)
#
#         # ====================
#         #   Noise constraint: ω ∈ {1, 2, 3, 4, 5}
#         if stage != 1
#             @noise(sp, Ω = 1:N, ω == Ω)
#         else
#             # First stage is wait-and-see
#             @constraint(sp, ω == 0)
#         end
#
#         # ====================
#         #   Dynamics: ensures x ∈ {1, 2, 3, 4, 5}
#         @constraint(sp, x == x0 + ω)
#
#         # ====================
#         #   Objective: maximise x, so future costs = [1 2 3 4 5]
#         stageobjective!(sp, x)
#
#     end
#
#     status = solve(m, max_iterations = 10)
#
#     return m, status
# end
#
# # ==================================
# # Build model with DRO, radius 0.2
# measure = DRO(0.2, eye(N))
# m, status = buildandsolvemodel(measure)
#
# @test isapprox(getbound(m), 2.367545, atol=1e-4)
#
# # ==================================
# # Build model with DRO, radius 0.25
# measure = DRO(0.25)
# m, status = buildandsolvemodel(measure)
#
# @test isapprox(getbound(m), 2.2094305, atol=1e-5)
#
# # ==================================
# # Build model with DRO, radius 0.4
# measure = DRO(0.4)
# m, status = buildandsolvemodel(measure)
#
# @test isapprox(getbound(m), 1.75838, atol=1e-5)
