using Kokako, Test, JSON, Gurobi, Plots

"""
    infinite_powder(; discount_factor = 0.5, stocking_rate::Float64 = NaN)

Create an instance of the infinite horizon POWDER model.
"""
function infinite_powder(; discount_factor = 0.5, stocking_rate::Float64 = NaN)
    data = JSON.parsefile(joinpath(@__DIR__, "powder_data.json"))
    # Allow over-ride of the stocking rate contained in data.
    if !isnan(stocking_rate)
        data["stocking_rate"] = stocking_rate
    end
    # ===== Markovian Graph =====
    transition = Array{Float64, 2}[]
    for transition_matrix in data["transition"]
        push!(
            transition,
            convert(
                Array{Float64, 2},
                reshape(
                    vcat(transition_matrix...),
                    length(transition_matrix[1]),
                    length(transition_matrix)
                )
            )
        )
    end
    graph = Kokako.MarkovianGraph(transition)
    for markov_state in 1:size(transition[end], 2)
        Kokako.add_edge(graph,
            (data["number_of_weeks"], markov_state) => (1, 1),
            discount_factor
        )
    end

    model = Kokako.PolicyGraph(graph,
        sense = :Max,
        bellman_function = Kokako.AverageCut(upper_bound = 1e5),
        optimizer = with_optimizer(Gurobi.Optimizer, OutputFlag = 0)
            ) do subproblem, index
        # Unpack the node index.
        stage, markov_state = index
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
        @variables(subproblem, begin
            # Pasture cover (kgDM/ha). Note: to avoid numerical difficulties, we
            # increase the lower bound so that it is not zero. This avoids the
            # situaton where pasture_cover=0 and thus growth=0, effectively
            # killing all grass for all time.
            (10 <= pasture_cover <= data["maximum_pasture_cover"], Kokako.State,
                initial_value = data["initial_pasture_cover"])
            # Quantity of supplement in storage (kgDM/ha).
            (stored_supplement >= 0, Kokako.State,
                initial_value = data["initial_storage"])
            # Soil moisture (mm).
            (0 <= soil_moisture <= data["maximum_soil_moisture"], Kokako.State,
                initial_value = data["initial_soil_moisture"])
            # Number of cows milking (cows/ha).
            (0 <= cows_milking <= data["stocking_rate"], Kokako.State,
                initial_value = data["stocking_rate"])
            (0 <= milk_production <= data["maximum_milk_production"],
                Kokako.State, initial_value = 0.0)
        end)
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
            weekly_milk_production >= 0  # Weekly milk production (kgMS/week).
            fei_penalty >= 0  # Fat Evaluation Index penalty ($)
        end)

        # ========== Parameterize model on uncertainty ==========
        Kokako.parameterize(subproblem, data["niwa_data"][stage]) do ω
            JuMP.set_upper_bound(evapotranspiration, ω["evapotranspiration"])
            JuMP.fix(rainfall, ω["rainfall"])
        end

        @constraints(subproblem, begin
            # ========== State constraints ==========
            pasture_cover.out ==
                pasture_cover.in + 7 * grass_growth - harvest - feed_pasture
            stored_supplement.out == stored_supplement.in +
                data["harvesting_efficiency"] * harvest - feed_storage
            # This is a <= do account for the maximum soil moisture; excess
            # water is assumed to drain away.
            soil_moisture.out <=
                soil_moisture.in - evapotranspiration + rainfall

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
            # Upper bound on the energy that can be used for milk production.
            energy_for_milk_production <=
                data["max_milk_energy"][stage] * cows_milking.in
            # Conversion between energy and physical milk
            weekly_milk_production == energy_for_milk_production /
                data["energy_content_of_milk"][stage]
            # Lower bound on milk production.
            weekly_milk_production >=
                data["min_milk_production"] * cows_milking.in

            # ========== Pasture growth models ==========
            # Model One: grass_growth ~ evapotranspiration
            grass_growth <=
                data["soil_fertility"][stage] * evapotranspiration / 7
            # Model Two: grass_growth ~ pasture_cover
            [p′ = range(0, stop = Pₘ, length = Pₙ)],
                grass_growth <= g(p′) + g′(p′) * (pasture_cover.in - p′)

            # ========== Fat Evaluation Index Penalty ==========
            fei_penalty >=
                cow_per_day * (0.00 + 0.25 * (supplement / cow_per_day - 3))
            fei_penalty >=
                cow_per_day * (0.25 + 0.50 * (supplement / cow_per_day - 4))
            fei_penalty >=
                cow_per_day * (0.75 + 1.00 * (supplement / cow_per_day - 5))
        end)

        # ========== Lactation cycle over the season ==========
        if stage == data["number_of_weeks"]
            @constraint(subproblem, cows_milking.out == data["stocking_rate"])
        elseif data["maximum_lactation"] <= stage < data["number_of_weeks"]
            @constraint(subproblem, cows_milking.out == 0)
        else
            @constraint(subproblem, cows_milking.out <= cows_milking.in)
        end

        # ========== Milk revenue cover penalty ==========
        if stage == data["number_of_weeks"]
            @constraint(subproblem, milk_production.out == 0.0)
            @expression(subproblem, milk_revenue,
                data["prices"][stage][markov_state] * milk_production.in)
        else
            @constraint(subproblem, milk_production.out ==
                milk_production.in + weekly_milk_production)
            @expression(subproblem, milk_revenue, 0.0)
        end

        # ========== Stage Objective ==========
        @stageobjective(subproblem,
            milk_revenue -
            data["supplement_price"] * supplement -
            data["harvest_cost"] * harvest -
            fei_penalty +
            # Artificial term to encourage max soil moisture.
            1e-4 * soil_moisture.out
        )
    end
    return model
end

function visualize_policy(model, filename)
    simulations = Kokako.simulate(model, 5_000, [
        :cows_milking,
        :pasture_cover,
        :soil_moisture,
        :grass_growth,
        :supplement,
        :weekly_milk_production,
        :fei_penalty
        ],
        terminate_on_cycle = false,
        max_depth = 52 * 5
    )
    open(filename * ".json", "w") do io
        write(io, JSON.json(simulations))
    end
    xticks = (1:26:5*52, repeat(["Aug", "Feb"], outer=5))
    plot(
        Kokako.publicationplot(simulations,
            data -> data[:cows_milking].out,
            title = "(a)",
            ylabel = "Cows Milking (cows/ha)",
            xticks = xticks),
        Kokako.publicationplot(simulations,
            data -> data[:pasture_cover].out / 1000,
            ylabel = "Pasture Cover (t/ha)",
            title = "(b)",
            xticks = xticks),
        Kokako.publicationplot(simulations,
            data -> data[:noise_term]["evapotranspiration"],
            ylabel = "Evapotranspiration (mm)",
            title = "(c)",
            xticks = xticks),
        # Kokako.publicationplot(simulations,
        #     data -> data[:grass_growth],
        #     ylabel = "Grass Growth (kg/day)",
        #     title = "(d)",
        #     xticks = xticks),
        # Kokako.publicationplot(simulations,
        #     data -> data[:supplement],
        #     ylabel = "Palm Kernel Fed (kg/cow/day)",
        #     title = "(e)",
        #     xticks = xticks),
        # Kokako.publicationplot(simulations,
        #     data -> data[:weekly_milk_production],
        #     ylabel = "Milk Production (kg/day)",
        #     title = "(f)",
        #     xticks = xticks),

        layout = (1, 3),
        size = (1500, 300)
    )
    savefig(filename * ".pdf")
end


# The experiments can be run by calling `julia powder.jl run`.
if length(ARGS) > 0 && ARGS[1] == "run"
    # Import Random for seed!
    using Random
    let
        println("""Running the base case: POWDER with discount factor of 0.9 and
        a stocking rate of 3 cows/ha.""")
        model = infinite_powder(discount_factor = 0.9, stocking_rate = 3.0)
        Random.seed!(1234)
        Kokako.train(model, iteration_limit = 3_000, print_level = 1,
            log_file = "powder_09_3000it.log")
        Random.seed!(5678)
        visualize_policy(model, "powder_09_3000it")
    end

    let
        println("""Running experiment 1: compare the evolution of the lower
        bound with a smaller # discount factor. No need to visualize as we only
        care about the evolution of the upper bound.""")
        model = infinite_powder(discount_factor = 0.5, stocking_rate = 3.0)
        Random.seed!(1234)
        Kokako.train(model, iteration_limit = 3_000, print_level = 1,
            log_file = "powder_50_3000it.log")
    end

    let
        println("""Running experiment 2: visualize an unconverged policy.
        Compare against powder_90_3000it.""")
        model = infinite_powder(discount_factor = 0.9, stocking_rate = 3.0)
        Random.seed!(1234)
        Kokako.train(model, iteration_limit = 500, print_level = 1,
            log_file = "powder_90_500it.log")
        Random.seed!(5678)
        visualize_policy(model, "powder_90_500it")
    end

    let
        println("""Running experiment 3: solve with a lower stocking rate.
        Compare against powder_90_3000it. No need to visualize as we only care
        about the lower bound.""")
        model = infinite_powder(discount_factor = 0.9, stocking_rate = 2.5)
        Random.seed!(1234)
        Kokako.train(model, iteration_limit = 3_000, print_level = 1,
            log_file = "powder_90_3000it_25SR.log")
    end
end
