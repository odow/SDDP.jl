using SDDP, JuMP, Clp
m = SDDPModel(
                 sense = :Min,
                stages = 3,
                solver = ClpSolver(),
       objective_bound = 0.0,
           is_infinite = true,
             lb_states = [0],
             ub_states = [200]) do sp, t

    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
    )
    setnoiseprobability!(sp, [1/3, 1/3, 1/3])
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end


function ff(x)
    s = SDDP.humanize(sum(x))
    while length(s) < 7
        s = string(" ", s)
    end
    s
end

solve(m, iteration_limit = 10, update_limit = 10, init_state_info_func=ff)
