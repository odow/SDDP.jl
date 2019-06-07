using SDDP, GLPK, Test

#=
TODO s

- user is reponsible for making all state variables binary

=#


# Number of time periods, number of seats, R_ij = evenue from selling seat i at time j, offer_ij = whether an offer for seat i will come at time j
(T, N, R, offer) = (3, 2, [3 3 6; 3 3 6], [1 1 0; 1 0 1])

model = SDDP.LinearPolicyGraph(
    stages = T,
    sense = :Max,
    upper_bound = 100.0,
    optimizer = with_optimizer(GLPK.Optimizer)) do sp, stage

    # Seat remaining?
    @variable(sp, 0 <= x[i in 1:N] <= 1, SDDP.State, initial_value = 1) # TODO why were bounds needed if binary?
    # Action: accept offer, or don't accept offer
    @variable(sp, accept_offer, Bin)
    # @variable(sp, 0 <= accept_offer <= 1)
    # Balance on seats
    @constraint(sp, [i in 1:N], x[i].out == x[i].in - offer[i, stage] * accept_offer)
    @stageobjective(sp, sum(R[i, stage] * offer[i, stage] * accept_offer for i in 1:N))
end

# n = model.nodes[2]
# SDDP.set_objective(n)
# sp = n.subproblem
# fix(sp[:x][1].in, 0.0)
# fix(sp[:x][2].in, 0.0)
# dual_vars = zeros(2)
# (obj, duals) = SDDP._kelley(n, dual_vars)
#
# n = model.nodes[3]
# SDDP.set_objective(n)
# sp = n.subproblem
# fix(sp[:x][1].in, 0.0)
# fix(sp[:x][2].in, 0.0)
# dual_vars = zeros(2)
# (obj, duals) = SDDP._kelley(n, dual_vars)


SDDP.train(model, iteration_limit = 20, print_level = 1)
@test SDDP.calculate_bound(model) ≈ 9.0


;
