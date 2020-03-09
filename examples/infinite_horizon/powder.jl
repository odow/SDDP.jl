#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, Test, JSON, Gurobi, Plots, Random

"""
    infinite_powder(; discount_factor = 0.75, stocking_rate::Float64 = NaN,
                    data_filename = "powder_data.json")

Create an instance of the infinite horizon POWDER model. If `stocking_rate =
NaN`, we use the value from the file `data_filename`.
"""
function infinite_powder(;
    discount_factor = 0.95,
    stocking_rate::Float64 = NaN,
    data_filename = "powder_data.json",
)
    data = JSON.parsefile(joinpath(@__DIR__, data_filename))
    # Allow over-ride of the stocking rate contained in data.
    if !isnan(stocking_rate)
        data["stocking_rate"] = stocking_rate
    end
    # ===== Markovian Graph =====
    transition = Array{Float64,2}[]
    for transition_matrix in data["transition"]
        push!(
            transition,
            convert(
                Array{Float64,2},
                Base.reshape(
                    vcat(transition_matrix...),
                    length(transition_matrix[1]),
                    length(transition_matrix),
                ),
            ),
        )
    end
    graph = SDDP.MarkovianGraph(transition)
    for markov_state = 1:size(transition[end], 2)
        SDDP.add_edge(
            graph,
            (data["number_of_weeks"], markov_state) => (1, 1),
            discount_factor,
        )
    end

    gurobi_env = Gurobi.Env()
    model = SDDP.PolicyGraph(
        graph,
        sense = :Max,
        upper_bound = 1e6,
        optimizer = () -> Gurobi.Optimizer(gurobi_env),
    ) do subproblem, index
        set_optimizer_attribute(subproblem, "OutputFlag", 0)
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
        @variables(
            subproblem,
            begin
                # Pasture cover (kgDM/ha). Note: to avoid numerical difficulties, we
                # increase the lower bound so that it is not zero. This avoids the
                # situaton where pasture_cover=0 and thus growth=0, effectively
                # killing all grass for all time.
                (
                    10 <= pasture_cover <= data["maximum_pasture_cover"],
                    SDDP.State,
                    initial_value = data["initial_pasture_cover"],
                )
                # Quantity of supplement in storage (kgDM/ha).
                (
                    stored_supplement >= 0,
                    SDDP.State,
                    initial_value = data["initial_storage"],
                )
                # Soil moisture (mm).
                (
                    0 <= soil_moisture <= data["maximum_soil_moisture"],
                    SDDP.State,
                    initial_value = data["initial_soil_moisture"],
                )
                # Number of cows milking (cows/ha).
                (
                    0 <= cows_milking <= data["stocking_rate"],
                    SDDP.State,
                    initial_value = data["stocking_rate"],
                )
                (
                    0 <= milk_production <= data["maximum_milk_production"],
                    SDDP.State,
                    initial_value = 0.0,
                )
            end
        )
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
        SDDP.parameterize(subproblem, data["niwa_data"][stage]) do ω
            JuMP.set_upper_bound(evapotranspiration, ω["evapotranspiration"])
            JuMP.fix(rainfall, ω["rainfall"])
        end

        @constraints(
            subproblem,
            begin
                # ========== State constraints ==========
                pasture_cover.out ==
                pasture_cover.in + 7 * grass_growth - harvest - feed_pasture
                stored_supplement.out ==
                stored_supplement.in + data["harvesting_efficiency"] * harvest -
                feed_storage
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
                # Upper bound on the energy that can be used for milk production.
                energy_for_milk_production <=
                data["max_milk_energy"][stage] * cows_milking.in
                # Conversion between energy and physical milk
                weekly_milk_production ==
                energy_for_milk_production / data["energy_content_of_milk"][stage]
                # Lower bound on milk production.
                weekly_milk_production >= data["min_milk_production"] * cows_milking.in

                # ========== Pasture growth models ==========
                # Model One: grass_growth ~ evapotranspiration
                grass_growth <= data["soil_fertility"][stage] * evapotranspiration / 7
                # Model Two: grass_growth ~ pasture_cover
                [p′ = range(0, stop = Pₘ, length = Pₙ)],
                grass_growth <= g(p′) + g′(p′) * (pasture_cover.in - p′)

                # ========== Fat Evaluation Index Penalty ==========
                fei_penalty >= cow_per_day * (0.00 + 0.25 * (supplement / cow_per_day - 3))
                fei_penalty >= cow_per_day * (0.25 + 0.50 * (supplement / cow_per_day - 4))
                fei_penalty >= cow_per_day * (0.75 + 1.00 * (supplement / cow_per_day - 5))
            end
        )

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
            @expression(
                subproblem,
                milk_revenue,
                data["prices"][stage][markov_state] * milk_production.in
            )
        else
            @constraint(
                subproblem,
                milk_production.out == milk_production.in + weekly_milk_production
            )
            @expression(subproblem, milk_revenue, 0.0)
        end

        # ========== Stage Objective ==========
        @stageobjective(
            subproblem,
            milk_revenue - data["supplement_price"] * supplement -
            data["harvest_cost"] * harvest - fei_penalty +
            # Artificial term to encourage max soil moisture.
            1e-4 * soil_moisture.out
        )
    end
    return model
end

function visualize_policy(model, filename)
    simulations = SDDP.simulate(
        model,
        1_000,
        [
            :cows_milking,
            :pasture_cover,
            :soil_moisture,
            :grass_growth,
            :supplement,
            :weekly_milk_production,
            :fei_penalty,
        ],
        sampling_scheme = SDDP.InSampleMonteCarlo(
            terminate_on_cycle = false,
            terminate_on_dummy_leaf = false,
            max_depth = 52 * 5,
        ),
    )
    xticks = (1:26:5*52, repeat(["Aug", "Feb"], outer = 5))
    plot(
        SDDP.publicationplot(
            simulations,
            data -> data[:cows_milking].out,
            title = "(a)",
            ylabel = "Cows Milking (cows/ha)",
            xticks = xticks,
        ),
        SDDP.publicationplot(
            simulations,
            data -> data[:pasture_cover].out / 1000,
            ylabel = "Pasture Cover (t/ha)",
            title = "(b)",
            xticks = xticks,
        ),
        SDDP.publicationplot(
            simulations,
            data -> data[:noise_term]["evapotranspiration"],
            ylabel = "Evapotranspiration (mm)",
            xlabel = " ",
            title = "(c)",
            xticks = xticks,
        ),
        layout = (1, 3),
        size = (1500, 300),
    )
    savefig(filename * ".pdf")
end

function estimate_statistical_bound(model, filename)
    # Simulate to estimate the lower (statistical) bound. Note that we need to
    # set `terminate_on_dummy_leaf = true`.
    bound_simulations = SDDP.simulate(
        model,
        1_000,
        sampling_scheme = SDDP.InSampleMonteCarlo(
            terminate_on_cycle = false,
            terminate_on_dummy_leaf = true,
        ),
    )
    objectives = [sum(x[:stage_objective] for x in sim) for sim in bound_simulations]

    open(filename * ".json", "w") do io
        write(io, JSON.json(objectives))
    end
end

# The experiments can be run by calling `julia powder.jl run`.
if length(ARGS) > 0
    if ARGS[1] == "run"
        model = infinite_powder(discount_factor = 0.95, stocking_rate = 3.0)
        Random.seed!(123)
        SDDP.train(
            model,
            iteration_limit = 1_000,
            print_level = 1,
            log_file = "powder_complete.log",
        )
        Random.seed!(456)
        visualize_policy(model, "powder_visualization")

        model = infinite_powder(discount_factor = 0.95, stocking_rate = 3.0)
        for loop = 1:5
            Random.seed!(123 * loop)
            SDDP.train(
                model,
                iteration_limit = 200,
                print_level = 1,
                log_file = "powder_$(loop).log",
            )
            Random.seed!(456 * loop)
            estimate_statistical_bound(model, "powder_bound_$(loop)")
        end
    elseif ARGS[1] == "summarize"
        using Statistics
        function modified_cox(X, α = 1.96)
            N = length(X)
            logX = log.(X)
            μ = Statistics.mean(logX)
            σ² = Statistics.var(logX)
            half_width = α * sqrt(σ² / N + σ²^2 / (2N - 2))
            return exp(μ + σ² / 2 - half_width), exp(μ + σ² / 2 + half_width)
        end
        function normal(X, α = 1.96)
            N = length(X)
            μ = Statistics.mean(X)
            σ = Statistics.std(X)
            return μ + α * σ / sqrt(N), μ - α * σ / sqrt(N)
        end
        for i = 1:5
            data = JSON.parsefile("powder_bound_$(i).json", use_mmap = false)
            println(i, " ", modified_cox(data))
        end
        for i = 1:5
            data = JSON.parsefile("powder_bound_$(i).json", use_mmap = false)
            println(i, " ", normal(data))
        end
    end
end
