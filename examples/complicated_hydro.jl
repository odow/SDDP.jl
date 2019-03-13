using Kokako, JSON, GLPK

const DATA = JSON.parsefile(joinpath(@__DIR__, "complicated_hydro.json"))

const T = 12
const PRICES = [18 + round(5 * sin(0.5 * (t - 2 - 1)), digits=2) for t in 1:T]
const FLOW_KNOTS = [50.0, 60.0, 70.0]
const POWER_KNOTS = [55.0, 65.0, 70.0]

model = Kokako.LinearPolicyGraph(
            stages = T,
            sense = :Min,
            bellman_function = Kokako.AverageCut(
                lower_bound = 0,
                deletion_minimum = 1_000_000
            ),
            optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, t
    @variable(subproblem, 0 <= volume[1:3] <= 200, Kokako.State, initial_value = 50)
    @variable(subproblem, inflow[i=1:3], Kokako.State, initial_value = DATA["initial_inflow"][1][i])
    @variables(subproblem, begin
        thermal_generation >= 0
        thermal_cost >= 0
        hydro_flow[1:3] >= 0
        hydro_spill[1:3] >= 0
        pour[1:3] >= 0
        hydro_generation >= 0
        0 <= dispatch[1:3, 1:3] <= 1
        ω[1:3]
    end)
    @constraints(subproblem, begin
        thermal_cost >= 10 * thermal_generation + 0
        thermal_cost >= 20 * thermal_generation - 500
        thermal_cost >= 50 * thermal_generation - 3_500
        volume[1].out == volume[1].in + inflow[1].out - hydro_flow[1] -
            hydro_spill[1] + pour[1]
        [i=2:3], volume[i].out == volume[i].in + inflow[i].out -
            hydro_flow[i] - hydro_spill[i] + pour[i] +
            hydro_flow[i-1] + hydro_spill[i-1]
        hydro_generation == sum(
            sum(POWER_KNOTS[j] * dispatch[i, j] for j in 1:3)
            for i in 1:3)
        [i=1:3], hydro_flow[i] == sum(FLOW_KNOTS[j] * dispatch[i, j] for j in 1:3)
        [i=1:3], sum(dispatch[i, j] for j in 1:3) <= 1
        hydro_generation + thermal_generation >= 600
    end)
    @stageobjective(subproblem,
        thermal_cost -
        PRICES[t] * (hydro_generation + thermal_generation) +
        10_000 * sum(pour)
    )
    if t == 1
        @constraint(subproblem, [i=1:3], volume[i].out >= 30)
        for i in 1:3
            JuMP.fix(inflow[i].out, DATA["initial_inflow"][2][i])
        end
    else
        for i in 1:3
            R = DATA["ar_matrix"]["$(t-1)"]["$(i-1)"]
            @constraint(subproblem,
                inflow[i].out ==
                    sum(get(R, "$(j-1)", 0.0) * inflow[j].in for j in 1:3) + ω[i])
        end
        Kokako.parameterize(subproblem, [1, 2]) do ϕ
            for i in 1:3
                JuMP.fix(ω[i], DATA["RHS_noise"][i][ϕ][t])
            end
        end
    end
    if t == T
        @constraint(subproblem, [i=1:3], volume[i].out >= 30)
    end
end

Kokako.train(model, iteration_limit = 50, print_level = 0)

@test Kokako.calculate_bound(model) ≈ 129_469 atol=1
