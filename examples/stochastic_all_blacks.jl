using SDDP, GLPK, Test
using Random
Random.seed!(11111)

T = 3 # Number of time periods
N = 2 # Number of seats
R = [3 3 6; 3 3 6] # R_ij = price of seat i at time j
s = 3 # Number of noises
offers = [[rand([0 1], N) for _ in 1:s] for _ in 1:T]

model = SDDP.LinearPolicyGraph(
    stages = T,
    sense = :Max,
    upper_bound = 100.0,
    optimizer = with_optimizer(GLPK.Optimizer),
    mip_solver = SDDP.SDDiP()) do sp, stage

    # Seat remaining?
    @variable(sp, 0 <= x[1:N] <= 1, SDDP.State, Bin, initial_value = 1) # TODO why were bounds needed if binary?
    # Action: accept offer, or don't accept offer
    # We are allowed to accpect some of the seats offered but not others in this formulation
    @variable(sp, accept_offer[1:N], Bin)
    @variable(sp, offers_made[1:N])
    # Balance on seats
    @constraint(sp, balance[i in 1:N], x[i].in - x[i].out == accept_offer[i])
    @stageobjective(sp, sum(R[i, stage] * accept_offer[i] for i in 1:N))
    SDDP.parameterize(sp, offers[stage]) do o
        JuMP.fix.(offers_made, o)
    end
    @constraint(sp, accept_offer .<= offers_made)

end

SDDP.train(model, iteration_limit = 10, print_level = 1)
@test SDDP.calculate_bound(model) ≈ 8.0
