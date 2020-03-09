#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This should kind of work, but it doesn't.
using SDDP, Ipopt, Test

function infinite_ball_on_beam()
    graph = SDDP.Graph(
        :root_node,
        [:time_step],
        [(:root_node => :time_step, 1.0), (:time_step => :time_step, 0.999)],
    )
    model = SDDP.PolicyGraph(
        graph,
        bellman_function = SDDP.BellmanFunction(lower_bound = 0.0),
        optimizer = Ipopt.Optimizer,
        direct_mode = false,
        sense = :Min,
    ) do subproblem, node
        set_optimizer_attribute(subproblem, "print_level", 0)
        Δt = 0.1  # time-step (s)
        m = 0.1  # mass (kg)
        J = 0.5  # moment of inertia (kg m²)
        g = 9.81  # graviational acceleration (m s⁻²)
        τ = 1  #3  # Maximum torque.
        @variables(subproblem, begin
            # Beam position.
            r, SDDP.State, (initial_value = 1)
            # Rate of change of beam position.
            r′, SDDP.State, (initial_value = 0)
            # Angle of beam.
            θ, SDDP.State, (initial_value = -0.1745)
            # Rate of change of angle of beam.
            θ′, SDDP.State, (initial_value = 0)

            # Control variable: torque to apply to beam.
            -τ <= u <= τ
        end)
        @constraints(subproblem, begin
            r.out - r.in == Δt * r′.in
            θ.out - θ.in == Δt * θ′.in
        end)
        @NLconstraints(
            subproblem,
            begin
                r′.out - r′.in == Δt * (r.in * θ′.in^2 - g * sin(θ.in))
                θ′.out - θ′.in ==
                Δt * (-2 * m * r.in * r′.in + -1 * m * g * r.in * cos(θ.in) + u) /
                (m * r.in^2 + J)
            end
        )
        @stageobjective(
            subproblem,
            100 * r.out^2 + r′.out^2 + θ.out^2 + θ′.out^2 + 0.01 * u^2
        )
    end
    return model
end

using Random
Random.seed!(1234)
model = infinite_ball_on_beam()

SDDP.train(
    model,
    iteration_limit = 2_000,
    print_level = 1,
    cycle_discretization_delta = 0.001,
)

maximum_depth = 0
for depth = 10:10:100
    try
        simulation = SDDP.simulate(model; max_depth = depth, terminate_on_cycle = false)
        global maximum_depth = depth
    catch
        print("Max_depth = $(maximum_depth)")
        break
    end
end

simulation = SDDP.simulate(
    model,
    1,
    [:r, :θ, :u],
    max_depth = maximum_depth,
    terminate_on_cycle = false,
)
using Plots
plot(
    plot(
        [s[:r].out for s in simulation[1]],
        ylabel = "Displacement",
        ylims = (-1, 3),
        xlims = (1, maximum_depth),
    ),
    plot(
        [s[:θ].out for s in simulation[1]],
        ylabel = "Angle",
        ylims = (-1, 1),
        xlims = (1, maximum_depth),
    ),
    plot(
        [s[:u] for s in simulation[1]],
        ylabel = "Torque",
        ylims = (-3, 3),
        xlims = (1, maximum_depth),
    ),
    layout = (3, 1),
    size = (1500, 500),
)
