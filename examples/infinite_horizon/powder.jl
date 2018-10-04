using Kokako, Test, JSON, Gurobi, Plots

function infinite_powder()
    data = JSON.parsefile(joinpath(@__DIR__, "powder_data.json"))
    # Construct infinite horizon graph.
    graph = Kokako.LinearGraph(data["number_of_weeks"])
    Kokako.add_edge(graph, data["number_of_weeks"] => 1, 0.5)

    model = Kokako.PolicyGraph(graph,
        sense = :Max,
        bellman_function = Kokako.AverageCut(upper_bound = 1e5),
        optimizer = with_optimizer(Gurobi.Optimizer, OutputFlag=0)
            ) do subproblem, stage
        # Simple linspace since it got dropped from Base Julia v1.0.
        linspace(a::Real, b::Real, N::Int) = collect(a:(b-a)/(N-1):b)
        # ========== Data Initialization ==========
        # Data for Fat Evaluation Index penalty
        cow_per_day = data["stocking_rate"] * 7
        # Data for grass growth model two
        Pₘ = data["maximum_pasture_cover"]
        gₘ = data["maximum_growth_rate"]
        Pₙ = data["number_of_pasture_cuts"]
        g(p) = 4 * gₘ / Pₘ * p * (1 - p / Pₘ)
        g′(p) = 4 * gₘ / Pₘ * (1 - 2 * p / Pₘ)

        # ========== State Variables ==========
        # Pasture cover (kgDM/ha).
        @variable(subproblem,
            0 <= pasture_cover <= data["maximum_pasture_cover"],
            Kokako.State, root_value = data["initial_pasture_cover"])
        # Quantity of supplement in storage (kgDM/ha).
        @variable(subproblem, stored_supplement >= 0,
            Kokako.State, root_value = data["initial_storage"])
        # Soil moisture (mm).
        @variable(subproblem, 0 <= soil_moisture <= data["maximum_soil_moisture"],
            Kokako.State, root_value = data["initial_soil_moisture"])
        # Number of cows milking (cows/ha).
        @variable(subproblem, 0 <= cows_milking <= data["stocking_rate"],
            Kokako.State, root_value = data["stocking_rate"])

        # ========== Control Variables ==========
        @variables(subproblem, begin
            supplement >= 0  # Quantity of supplement to buy and feed (kgDM).
            harvest >= 0  # Quantity of pasture to harvest (kgDM/ha).
            feed_storage >= 0  # Feed herd grass from storage (kgDM).
            feed_pasture >= 0  # Feed herd grass from pasture (kgDM).
            evapotranspiration >= 0  # The actual evapotranspiration rate.
            rainfall  # Rainfall (mm); dummy variable for parameterization.
            grass_growth >= 0  # The potential grass growth rate.
            energy_for_milk_production >= 0  # Energy for milk production (MJ).
            actual_milk_production >= 0  # Actual milk production (kgMS).
            fei_penalty >= 0  # Fat Evaluation Index penalty ($)
        end)

        # ========== Parameterize model on uncertainty ==========
        Kokako.parameterize(subproblem, data["niwa_data"][stage]) do ω
            JuMP.set_upper_bound(evapotranspiration, ω["evapotranspiration"])
            JuMP.fix(rainfall, ω["rainfall"])
        end

        @constraints(subproblem, begin
            # ========== State constraints ==========
            pasture_cover.out <= pasture_cover.in + 7 * grass_growth - harvest - feed_pasture
            stored_supplement.out <= stored_supplement.in + data["harvesting_efficiency"] * harvest - feed_storage
            # This is a <= do account for the maximum soil moisture; excess
            # water is assumed to drain away.
            soil_moisture.out <= soil_moisture.in - evapotranspiration + rainfall

            # ========== Energy balance ==========
            data["pasture_energy_density"] * (feed_pasture + feed_storage) +
                data["supplement_energy_density"] * supplement >=
                data["stocking_rate"] * (
                    data["energy_for_pregnancy"][stage] +
                    data["energy_for_maintenance"] +
                    data["energy_for_bcs_dry"][stage]
                ) +
                cows_milking.in * (
                    data["energy_for_bcs_milking"][stage] -
                    data["energy_for_bcs_dry"][stage]
                ) +
                energy_for_milk_production

            # ========== Milk production models ==========
            energy_for_milk_production <= data["max_milk_energy"][stage] * cows_milking.in
            actual_milk_production <= energy_for_milk_production / data["energy_content_of_milk"][stage]
            # minimum milk
            actual_milk_production >= data["min_milk_production"] * cows_milking.in

            # ========== Pasture growth models ==========
            # Model One: grass_growth ~ evapotranspiration
            grass_growth <= data["soil_fertility"][stage] * evapotranspiration / 7
            # Model Two: grass_growth ~ pasture_cover
            [p′ = linspace(0, Pₘ, Pₙ)], grass_growth <= g(p′) + g′(p′) * (pasture_cover.in - p′)

            # ========== Fat Evaluation Index Penalty ==========
            fei_penalty >= cow_per_day * (0.00 + 0.25 * (supplement / cow_per_day - 3))
            fei_penalty >= cow_per_day * (0.25 + 0.50 * (supplement / cow_per_day - 4))
            fei_penalty >= cow_per_day * (0.75 + 1.00 * (supplement / cow_per_day - 5))
        end)

        # ========== Lactation cycle over the season ==========
        if stage == data["number_of_weeks"]
            @constraint(subproblem, cows_milking.out == data["stocking_rate"])
        elseif data["maximum_lactation"] <= stage < data["number_of_weeks"]
            @constraint(subproblem, cows_milking.out == 0)
        else
            @constraint(subproblem, cows_milking.out <= cows_milking.in)
        end

        # ========== Low pasture cover penalty ==========
        # To encourage the optimization to avoid zero pasture cover (thereby
        # enforcing grass growth to zero and all manner of numerical issues),
        # add a small penalty when the pasture over gets unreasonably low. We
        # should never see pasture_penalty > 0 in a simulation.
        @variable(subproblem, pasture_penalty >= 0)
        @constraint(subproblem, pasture_cover.out + pasture_penalty >= 500.0)

        # ========== Stage Objective ==========
        @stageobjective(subproblem,
            6.0 * actual_milk_production -
            pasture_penalty -
            data["supplement_price"] * supplement -
            data["harvest_cost"] * harvest -
            fei_penalty +
            1e-4 * soil_moisture.out  # Artificial term to encourage max soil moisture.
        )
    end

    Kokako.train(model, iteration_limit = 1_000, print_level = 1)

    simulations = Kokako.simulate(model, 500, [
        :cows_milking,
        :pasture_cover,
        :soil_moisture,
        :grass_growth,
        :supplement,
        :actual_milk_production,
        :fei_penalty
        ],
        terminate_on_cycle = false,
        max_depth = 520
    )

    return model, simulations
end

model, simulations = infinite_powder()

plot(
    Kokako.publicationplot(simulations,
        data -> data[:cows_milking].out,
        ylabel = "Cows Milking (cows/ha)"),
    Kokako.publicationplot(simulations,
        data -> data[:pasture_cover].out / 1000,
        ylabel = "Pasture Cover (t/ha)"),
    Kokako.publicationplot(simulations,
        data -> data[:soil_moisture].out,
        ylabel = "Soil Moisture (mm)"),

    Kokako.publicationplot(simulations,
        data -> data[:grass_growth],
        ylabel = "Grass Growth (kg/day)"),
    Kokako.publicationplot(simulations,
        data -> data[:supplement],
        ylabel = "Palm Kernel Fed (kg/cow/day)"),
    Kokako.publicationplot(simulations,
        data -> data[:actual_milk_production],
        ylabel = "Milk Production (kg/day)"),

    Kokako.publicationplot(simulations,
        data -> data[:stage_objective],
        ylabel = "Stage Objective (\\\$)"),
    Kokako.publicationplot(simulations,
        data -> data[:noise_term]["evapotranspiration"],
        ylabel = "Evapotranspiration (mm)"),
    Kokako.publicationplot(simulations,
        data -> data[:noise_term]["rainfall"],
        ylabel = "Rainfall (mm)"),
    size = (1500, 900)
)
