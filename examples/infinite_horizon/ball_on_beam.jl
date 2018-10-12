using Kokako, Ipopt, Test

function value(ex::JuMP.GenericQuadExpr{CoefType, VarType},
               map::Function) where {CoefType, VarType}
    # The return type of appling map to ::VarType.
    MapVarType = Base.promote_op(map, VarType)
    # Later, we're going to multiply two MapVarType together
    MapVarType2 = Base.promote_op(*, MapVarType, MapVarType)
    # We're also going to multiple a constant with ::MapVarType2
    RetType = Base.promote_op(*, CoefType, MapVarType2)
    ret = convert(RetType, map(ex.aff))
    for (vars, coef) in ex.terms
        ret += coef * map(vars.a) * map(vars.b)
    end
    return ret
end

JuMP.result_value(ex::JuMP.GenericQuadExpr) = value(ex, JuMP.result_value)

function infinite_ball_on_beam()
    graph = Kokako.Graph(
        :root_node,
        [:time_step],
        [(:root_node => :time_step, 1.0), (:time_step => :time_step, 0.999)]
    )
    model = Kokako.PolicyGraph(graph,
                bellman_function = Kokako.AverageCut(lower_bound = 0.0),
                optimizer = with_optimizer(Ipopt.Optimizer, print_level = 0),
                direct_mode = false,
                sense = :Min
                    ) do subproblem, node
        Δt = 0.1  # time-step (s)
        m = 0.1  # mass (kg)
        J = 0.5  # moment of inertia (kg m²)
        g = 9.81  # graviational acceleration (m s⁻²)
        τ = 1  #3  # Maximum torque.
        @variables(subproblem, begin
            # Beam position.
            r, Kokako.State, (initial_value = 1)
            # Rate of change of beam position.
            r′, Kokako.State, (initial_value = 0)
            # Angle of beam.
            θ, Kokako.State, (initial_value = -0.1745)
            # Rate of change of angle of beam.
            θ′, Kokako.State, (initial_value = 0)

            # Control variable: torque to apply to beam.
            -τ <= u <= τ
        end)
        @constraints(subproblem, begin
            r.out - r.in == Δt * r′.in
            θ.out - θ.in == Δt * θ′.in
        end)
        @NLconstraints(subproblem, begin
            r′.out - r′.in == Δt * (r.in * θ′.in^2 - g * sin(θ.in))
            θ′.out - θ′.in == Δt * (
                -2 * m * r.in * r′.in +
                -1 * m * g * r.in * cos(θ.in) +
                u)  / (m * r.in^2 + J)
        end)
        @stageobjective(subproblem,
            100 * r.out^2 + r′.out^2 + θ.out^2 + θ′.out^2 + 0.01 * u^2)
    end
    return model
end

using Random
Random.seed!(1234)
model = infinite_ball_on_beam()

Kokako.train(model, iteration_limit = 2_000, print_level = 1,
    cycle_discretization_delta = 0.001)

maximum_depth = 0
for depth in 10:10:100
    try
        simulation = Kokako.simulate(model;
            max_depth = depth, terminate_on_cycle = false)
        global maximum_depth = depth
    catch
        print("Max_depth = $(maximum_depth)")
        break
    end
end

simulation = Kokako.simulate(model, 1, [:r, :θ, :u],
    max_depth = maximum_depth, terminate_on_cycle = false)
using Plots
plot(
    plot([s[:r].out for s in simulation[1]],
        ylabel = "Displacement", ylims = (-1, 3), xlims = (1, maximum_depth)),
    plot([s[:θ].out for s in simulation[1]],
        ylabel = "Angle", ylims = (-1, 1), xlims = (1, maximum_depth)),
    plot([s[:u] for s in simulation[1]],
        ylabel = "Torque", ylims = (-3, 3), xlims = (1, maximum_depth)),
    layout = (3, 1),
    size = (1500, 500)
)
