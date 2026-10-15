# Copyright (c) 2016-26: Electric Power Optimization Centre, University of
# Auckland, and Oscar Dowson and contributors.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

module TestJADE

using Test

import DBInterface
import HiGHS
import SDDP
import SQLite

# Load the JADE module.
include(joinpath(@__DIR__, "JADE.jl"))

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

const CSV_DIR = joinpath(@__DIR__, "data")

# The tests read a database that is converted from the CSV files once, in a
# temporary directory.
const DB_FILE = joinpath(mktempdir(), "data.db")
JADE.Data(; filename = DB_FILE, directory = CSV_DIR)

function test_data()
    data = JADE.Data(DB_FILE)
    @test data isa JADE.Data
    @test fieldnames(JADE.Data) == Tuple(Symbol.(JADE.TABLES))
    @test sprint(show, data) == "JADE.Data with 15 tables"
    # There is no constructor that takes the tables directly
    @test length(methods(JADE.Data)) == 2
    @test_throws MethodError JADE.Data(data.hours_per_block)
    @test size(data.demand) == (3 * 52 * 3, 4)
    @test size(data.hours_per_block) == (52 * 3, 3)
    @test size(data.inflows) == (26 * 50 * 52, 4)
    @test size(data.transmission) == (4 * 5, 5)
    # NULL means there is no limit
    row = only(filter(r -> r.ORIG == "Lake_Wanaka", data.hydro_arcs))
    @test ismissing(row.MAX_FLOW)
    return
end

function test_data_missing_file()
    @test_throws AssertionError JADE.Data("not_a_file.db")
    return
end

function test_data_wrong_schema()
    mktempdir() do dir
        filename = joinpath(dir, "empty.db")
        db = SQLite.DB(filename)
        DBInterface.execute(db, "CREATE TABLE demand (NODE TEXT)")
        DBInterface.close!(db)
        @test_throws AssertionError JADE.Data(filename)
    end
    return
end

function test_schema_rejects_invalid_rows()
    db = SQLite.DB()
    JADE._create_tables(db)
    insert = "INSERT INTO station_outages (GENERATOR, WEEK, OUTAGE) VALUES"
    DBInterface.execute(db, "$insert ('Huntly', 1, 250)")
    # Duplicate primary key
    @test_throws(
        SQLite.SQLiteException,
        DBInterface.execute(db, "$insert ('Huntly', 1, 100)"),
    )
    # WEEK must be in 1:52
    @test_throws(
        SQLite.SQLiteException,
        DBInterface.execute(db, "$insert ('Huntly', 53, 250)"),
    )
    # OUTAGE must be non-negative
    @test_throws(
        SQLite.SQLiteException,
        DBInterface.execute(db, "$insert ('Huntly', 2, -1)"),
    )
    # Tables are STRICT, so OUTAGE must be a number
    @test_throws(
        SQLite.SQLiteException,
        DBInterface.execute(db, "$insert ('Huntly', 3, 'large')"),
    )
    # Every table has a primary key, so duplicate rows are rejected
    insert = "INSERT INTO lost_load (NODE, COST, PROPORTION) VALUES"
    DBInterface.execute(db, "$insert ('NI', 1000, 0.5)")
    @test_throws(
        SQLite.SQLiteException,
        DBInterface.execute(db, "$insert ('NI', 1000, 0.5)"),
    )
    DBInterface.close!(db)
    return
end

function test_data_from_csv()
    mktempdir() do dir
        filename = joinpath(dir, "default.db")
        data = JADE.Data(; filename, directory = CSV_DIR)
        @test isfile(filename)
        # Converting the CSV files is deterministic
        default = JADE.Data(DB_FILE)
        for table in fieldnames(JADE.Data)
            @test isequal(
                sort(getfield(data, table)),
                sort(getfield(default, table)),
            )
        end
        # Don't overwrite an existing file unless forced
        @test_throws(
            AssertionError,
            JADE.Data(; filename, directory = CSV_DIR),
        )
        mtime_before = mtime(filename)
        data = JADE.Data(; filename, directory = CSV_DIR, force = true)
        @test mtime(filename) >= mtime_before
        @test isequal(sort(data.demand), sort(default.demand))
        # Without a directory, the file is read, not written
        @test_throws(
            AssertionError,
            JADE.Data(; filename = joinpath(dir, "missing.db")),
        )
        @test JADE.Data(; filename) isa JADE.Data
    end
    return
end

function test_csv_fixed_stations()
    df = JADE._csv_fixed_stations(CSV_DIR)
    @test names(df) == ["STATION", "NODE"]
    @test size(df, 1) == 3
    @test only(filter(r -> r.STATION == "WestWind", df)).NODE == "HAY"
    return
end

function test_csv_fixed_station_output()
    df = JADE._csv_fixed_station_output(CSV_DIR)
    @test names(df) == ["STATION", "LOAD_BLOCK", "OUTPUT"]
    @test size(df, 1) == 3 * 3
    row =
        only(filter(r -> r.STATION == "WestWind" && r.LOAD_BLOCK == "PEAK", df))
    @test row.OUTPUT == 60
    return
end

function test_csv_fuels()
    df = JADE._csv_fuels(CSV_DIR)
    @test names(df) == ["FUEL", "CO2_CONTENT"]
    @test df.FUEL == ["coal", "diesel", "gas"]
    @test df.CO2_CONTENT == [0.0912, 0.0, 0.0528]
    return
end

function test_csv_fuel_prices()
    df = JADE._csv_fuel_prices(CSV_DIR)
    @test names(df) == ["WEEK", "FUEL", "PRICE"]
    @test size(df, 1) == 52 * 3
    row = only(filter(r -> r.WEEK == 1 && r.FUEL == "gas", df))
    @test row.PRICE == 6.68
    return
end

function test_csv_co2_prices()
    df = JADE._csv_co2_prices(CSV_DIR)
    @test names(df) == ["WEEK", "CO2_PRICE"]
    @test df.WEEK == 1:52
    @test all(df.CO2_PRICE .== 25.0)
    return
end

function test_csv_hydro_arcs()
    df = JADE._csv_hydro_arcs(CSV_DIR)
    row = only(filter(r -> r.ORIG == "Lake_Wanaka", df))
    @test row.MIN_FLOW == 0.0
    @test ismissing(row.MAX_FLOW)
    # Spillways are hydro arcs
    row = only(filter(r -> r.ORIG == "Lake_Aviemore", df))
    @test row.DEST == "Lake_Waitaki"
    @test row.MAX_FLOW == 5400
    return
end

function test_csv_inflows()
    df = JADE._csv_inflows(CSV_DIR)
    @test names(df) == ["CATCHMENT", "YEAR", "WEEK", "INFLOW"]
    @test size(df, 1) == 26 * 50 * 52
    rows = filter(r -> r.CATCHMENT == "Lake_Aviemore" && r.YEAR == 1970, df)
    @test rows.WEEK[1:2] == [1, 2]
    @test rows.INFLOW[1:2] == [27, 16]
    return
end

function test_csv_lost_load()
    df = JADE._csv_lost_load(CSV_DIR)
    @test names(df) == ["NODE", "COST", "PROPORTION"]
    @test sort(df.NODE) == ["HAY", "NI", "SI"]
    @test all(df.PROPORTION .== 1.0)
    @test all(df.COST .== 1000)
    return
end

function test_csv_reservoirs()
    df = JADE._csv_reservoirs(CSV_DIR)
    @test names(df) == ["RESERVOIR", "INITIAL_LEVEL", "MAX_LEVEL"]
    row = only(filter(r -> r.RESERVOIR == "Lake_Taupo", df))
    @test row.INITIAL_LEVEL == 750.275
    @test row.MAX_LEVEL == 848.62423
    return
end

function test_csv_station_outages()
    df = JADE._csv_station_outages(CSV_DIR)
    @test names(df) == ["GENERATOR", "WEEK", "OUTAGE"]
    @test sort(unique(df.GENERATOR)) == ["Huntly_main_g3", "Huntly_main_g4"]
    @test size(df, 1) == 2 * 52
    @test all(df.OUTAGE .== 250)
    return
end

function test_csv_transmission()
    df = JADE._csv_transmission(CSV_DIR)
    @test names(df) ==
          ["FROM_NODE", "TO_NODE", "TRANCHE", "CAPACITY", "LOSS_FRACTION"]
    rows = filter(r -> r.FROM_NODE == "NI" && r.TO_NODE == "HAY", df)
    @test rows.TRANCHE == 1:5
    @test rows.CAPACITY == [100.0, 100.0, 100.0, 100.0, 1943.0]
    @test rows.LOSS_FRACTION == [0.005, 0.015, 0.025, 0.035, 0.05]
    @test sum(rows.CAPACITY) == 2343
    return
end

function test_transmission_losses()
    data = JADE.Data(DB_FILE)
    @test JADE._transmission_lines(data) ==
          [("HAY", "NI"), ("HAY", "SI"), ("NI", "HAY"), ("SI", "HAY")]
    @test JADE._transmission_losses(data, ("NI", "HAY")) == [
        (100.0, 0.005),
        (100.0, 0.015),
        (100.0, 0.025),
        (100.0, 0.035),
        (1943.0, 0.05),
    ]
    return
end

# Write `data` to a new database, with the tables in `kwargs` replaced, and
# read it back. Reading it back runs `_validate_data`.
function _data_with_tables(data::JADE.Data; kwargs...)
    filename = joinpath(mktempdir(), "data.db")
    db = SQLite.DB(filename)
    JADE._create_tables(db)
    for table in JADE.TABLES
        df = get(kwargs, Symbol(table), getfield(data, Symbol(table)))
        SQLite.load!(df, db, table)
    end
    DBInterface.close!(db)
    return JADE.Data(filename)
end

function test_validate_data()
    data = JADE.Data(DB_FILE)
    @test _data_with_tables(data) isa JADE.Data
    lost_load = copy(data.lost_load)
    lost_load.PROPORTION .= 0.5
    @test_throws AssertionError _data_with_tables(data; lost_load)
    station_outages = copy(data.station_outages)
    station_outages.OUTAGE .= 1_000.0
    @test_throws AssertionError _data_with_tables(data; station_outages)
    hydro_arcs = filter(r -> r.ORIG != "Lake_Waikaremoana", data.hydro_arcs)
    @test_throws AssertionError _data_with_tables(data; hydro_arcs)
    fuel_prices = filter(r -> r.FUEL != "gas", data.fuel_prices)
    @test_throws AssertionError _data_with_tables(data; fuel_prices)
    fixed_station_output =
        filter(r -> r.STATION != "WestWind", data.fixed_station_output)
    @test_throws AssertionError _data_with_tables(data; fixed_station_output)
    transmission = copy(data.transmission)
    transmission.LOSS_FRACTION .= reverse(transmission.LOSS_FRACTION)
    @test_throws AssertionError _data_with_tables(data; transmission)
    transmission = filter(r -> r.TRANCHE != 2, data.transmission)
    @test_throws AssertionError _data_with_tables(data; transmission)
    inflows = filter(r -> r.YEAR != 1990, data.inflows)
    @test_throws AssertionError _data_with_tables(data; inflows)
    return
end

function test_max_specific_power()
    data = JADE.Data(DB_FILE)
    specific_power = JADE._specific_power(data)
    @test specific_power["SEA"] == 0.0
    @test all(v -> v >= 0, values(specific_power))
    # Lake_Cobb flows straight through the Cobb station to the sea.
    @test specific_power["Lake_Cobb"] == 4.405
    @test JADE._max_specific_power(data) == 4.405
    return
end

function test_config()
    config = JADE.Config(; optimizer = HiGHS.Optimizer)
    @test config.optimizer === HiGHS.Optimizer
    @test config.current_week == 1
    @test config.discount_factor == 0.95
    @test config.penalty_below_min_flow == 500.0
    @test config.penalty_above_max_flow == 50.0
    @test config.objective_scale == 1e-6
    @test config.inflow_years === nothing
    @test_throws UndefKeywordError JADE.Config()
    # Abstract vectors are converted to a `Vector{Int}`
    config = JADE.Config(; optimizer = HiGHS.Optimizer, inflow_years = 1:3)
    @test config.inflow_years isa Vector{Int}
    @test config.inflow_years == [1, 2, 3]
    years = [2000, 1990]
    config = JADE.Config(; optimizer = HiGHS.Optimizer, inflow_years = years)
    @test config.inflow_years === years
    return
end

function test_policy_graph()
    data = JADE.Data(DB_FILE)
    config = JADE.Config(; optimizer = HiGHS.Optimizer)
    model = SDDP.PolicyGraph(data, config)
    @test model isa SDDP.PolicyGraph
    @test sort(collect(keys(model.nodes))) == 1:52
    @test model.root_children == [SDDP.Noise(1, 1.0)]
    @test length(model[1].noise_terms) == 50
    @test model[52].children == [SDDP.Noise(1, 0.95)]
    # The factors that convert duals into prices are stored in each subproblem
    sp = model[23].subproblem
    hours = JADE._hours(data, 23)
    @test sp.ext[:price_scale]["PEAK"] ≈ 1 / (1e-6 * hours["PEAK"])
    specific_power = JADE._specific_power(data)
    r = first(data.reservoirs.RESERVOIR)
    @test sp.ext[:water_value_scale][r] ≈
          -JADE.MM3_PER_CUMEC_HOUR / (1e-6 * specific_power[r])
    SDDP.train(
        model;
        iteration_limit = 2,
        sampling_scheme = SDDP.InSampleMonteCarlo(; max_depth = 52),
        print_level = 0,
    )
    @test SDDP.calculate_bound(model) > 0
    simulations = SDDP.simulate(
        model,
        2,
        [:x_storage];
        sampling_scheme = SDDP.InSampleMonteCarlo(;
            max_depth = 52,
            terminate_on_dummy_leaf = false,
        ),
    )
    @test length(simulations) == 2
    @test all(length.(simulations) .== 52)
    @test [s[:node_index] for s in simulations[1]] == 1:52
    return
end

function test_policy_graph_inflow_years()
    data = JADE.Data(DB_FILE)
    # A year is in-sample if and only if it is in `inflow_years`. The order
    # and duplicates in `inflow_years` do not matter.
    for inflow_years in ([1995], 1980:1984, [2001, 1973, 2019, 1973])
        config = JADE.Config(; optimizer = HiGHS.Optimizer, inflow_years)
        model = SDDP.PolicyGraph(data, config)
        years = sort(unique(inflow_years))
        for week in (1, 23, 52)
            noises = model[week].noise_terms
            @test length(noises) == length(years)
            @test all(noise.probability ≈ 1 / length(years) for noise in noises)
            terms = [_historical_inflows(data, year, week) for year in years]
            @test [noise.term for noise in noises] == terms
        end
    end
    # Every year in `inflow_years` must be in the data.
    config = JADE.Config(; optimizer = HiGHS.Optimizer, inflow_years = Int[])
    @test_throws AssertionError SDDP.PolicyGraph(data, config)
    for inflow_years in ([1969], 2015:2020, [1990, 2050])
        config = JADE.Config(; optimizer = HiGHS.Optimizer, inflow_years)
        @test_throws AssertionError SDDP.PolicyGraph(data, config)
    end
    return
end

function test_policy_graph_current_week()
    data = JADE.Data(DB_FILE)
    config = JADE.Config(;
        optimizer = HiGHS.Optimizer,
        current_week = 23,
        discount_factor = 0.9,
    )
    model = SDDP.PolicyGraph(data, config)
    @test model.root_children == [SDDP.Noise(23, 1.0)]
    @test model[22].children == [SDDP.Noise(23, 0.9)]
    @test model[52].children == [SDDP.Noise(1, 1.0)]
    return
end

function _historical_inflows(data, year, week)
    df = filter(r -> r.YEAR == year && r.WEEK == week, data.inflows)
    return Dict(df.CATCHMENT .=> df.INFLOW)
end

function test_historical_years()
    data = JADE.Data(DB_FILE)
    config = JADE.Config(; optimizer = HiGHS.Optimizer, current_week = 39)
    sampling_scheme = JADE.HistoricalYears(data, config)
    @test sampling_scheme isa SDDP.Historical
    scenarios = [noise.term for noise in sampling_scheme.scenarios]
    # 1970 to 2018: the last year cannot start a scenario
    @test length(scenarios) == 49
    for scenario in scenarios
        @test first.(scenario) == [39:52; 1:38]
    end
    # The first scenario is week 39 of 1970 to week 38 of 1971
    @test scenarios[1][1][2] == _historical_inflows(data, 1970, 39)
    @test scenarios[1][14][2] == _historical_inflows(data, 1970, 52)
    @test scenarios[1][15][2] == _historical_inflows(data, 1971, 1)
    @test scenarios[1][end][2] == _historical_inflows(data, 1971, 38)
    # The last scenario is week 39 of 2018 to week 38 of 2019
    @test scenarios[end][1][2] == _historical_inflows(data, 2018, 39)
    @test scenarios[end][end][2] == _historical_inflows(data, 2019, 38)
    # If the model starts in week 1, every year is a scenario
    config_1 = JADE.Config(; optimizer = HiGHS.Optimizer)
    sampling_scheme = JADE.HistoricalYears(data, config_1)
    scenarios = [noise.term for noise in sampling_scheme.scenarios]
    @test length(scenarios) == 50
    for scenario in scenarios
        @test first.(scenario) == 1:52
    end
    @test scenarios[1][1][2] == _historical_inflows(data, 1970, 1)
    @test scenarios[end][end][2] == _historical_inflows(data, 2019, 52)
    # The scenarios can be simulated by the model
    model = SDDP.PolicyGraph(data, config)
    simulations = SDDP.simulate(model, 2, [:x_storage]; sampling_scheme)
    @test length(simulations) == 2
    @test all(s -> length(s) == 52, simulations)
    return
end

function test_historical_sequence()
    data = JADE.Data(DB_FILE)
    config = JADE.Config(; optimizer = HiGHS.Optimizer, current_week = 39)
    sampling_scheme = JADE.HistoricalSequence(data, config; start_year = 1980)
    @test sampling_scheme isa SDDP.Historical
    scenario = only(sampling_scheme.scenarios).term
    # Week 39 of 1980 to week 38 of 2019
    @test length(scenario) == 52 * (2019 - 1980)
    @test scenario[1] == (39, _historical_inflows(data, 1980, 39))
    @test scenario[14] == (52, _historical_inflows(data, 1980, 52))
    @test scenario[15] == (1, _historical_inflows(data, 1981, 1))
    @test scenario[end] == (38, _historical_inflows(data, 2019, 38))
    # The default start year is the first historical year
    sampling_scheme = JADE.HistoricalSequence(data, config)
    scenario_default = only(sampling_scheme.scenarios).term
    @test length(scenario_default) == 52 * (2019 - 1970)
    @test scenario_default[1] == (39, _historical_inflows(data, 1970, 39))
    @test scenario_default[end] == (38, _historical_inflows(data, 2019, 38))
    # If the model starts in week 1, the scenario ends in week 52
    config_1 = JADE.Config(; optimizer = HiGHS.Optimizer)
    sampling_scheme = JADE.HistoricalSequence(data, config_1; start_year = 2018)
    scenario_1 = only(sampling_scheme.scenarios).term
    @test length(scenario_1) == 52 * 2
    @test scenario_1[1] == (1, _historical_inflows(data, 2018, 1))
    @test scenario_1[end] == (52, _historical_inflows(data, 2019, 52))
    # `end_year` stops the scenario in week 38 of `end_year`
    sampling_scheme = JADE.HistoricalSequence(
        data,
        config;
        start_year = 1980,
        end_year = 1990,
    )
    scenario_end = only(sampling_scheme.scenarios).term
    @test length(scenario_end) == 52 * (1990 - 1980)
    @test scenario_end[1] == (39, _historical_inflows(data, 1980, 39))
    @test scenario_end[end] == (38, _historical_inflows(data, 1990, 38))
    # If the model starts in week 1, a single year is possible
    sampling_scheme = JADE.HistoricalSequence(
        data,
        config_1;
        start_year = 2018,
        end_year = 2018,
    )
    scenario_end_1 = only(sampling_scheme.scenarios).term
    @test length(scenario_end_1) == 52
    @test scenario_end_1[1] == (1, _historical_inflows(data, 2018, 1))
    @test scenario_end_1[end] == (52, _historical_inflows(data, 2018, 52))
    # Invalid start and end years
    @test_throws(
        AssertionError,
        JADE.HistoricalSequence(data, config; start_year = 1900),
    )
    @test_throws(
        AssertionError,
        JADE.HistoricalSequence(data, config; end_year = 2020),
    )
    # Week 39 of 2019 would need inflows from 2020
    @test_throws(
        AssertionError,
        JADE.HistoricalSequence(data, config; start_year = 2019),
    )
    # Week 39 of 1980 to week 38 of 1980 is not a whole year
    @test_throws(
        AssertionError,
        JADE.HistoricalSequence(
            data,
            config;
            start_year = 1980,
            end_year = 1980,
        ),
    )
    # `end_year` must not be before `start_year`
    @test_throws(
        AssertionError,
        JADE.HistoricalSequence(
            data,
            config_1;
            start_year = 1990,
            end_year = 1980,
        ),
    )
    # The scenario can be simulated by the model
    model = SDDP.PolicyGraph(data, config)
    sampling_scheme = JADE.HistoricalSequence(data, config; start_year = 1980)
    simulations = SDDP.simulate(model, 1, [:x_storage]; sampling_scheme)
    @test length(only(simulations)) == length(scenario)
    @test [s[:node_index] for s in only(simulations)] == first.(scenario)
    @test [s[:noise_term] for s in only(simulations)] == last.(scenario)
    return
end

function test_simulate()
    data = JADE.Data(DB_FILE)
    config = JADE.Config(; optimizer = HiGHS.Optimizer)
    model = SDDP.PolicyGraph(data, config)
    SDDP.train(
        model;
        iteration_limit = 2,
        sampling_scheme = SDDP.InSampleMonteCarlo(; max_depth = 52),
        print_level = 0,
    )
    sampling_scheme = JADE.HistoricalSequence(data, config; start_year = 2019)
    simulations = JADE.simulate(model, 1; sampling_scheme)
    simulation = only(simulations)
    @test length(simulation) == 52
    @test [s[:node_index] for s in simulation] == 1:52
    R, G = data.reservoirs.RESERVOIR, data.thermal_stations.GENERATOR
    P, B = unique(data.demand.NODE), unique(data.hours_per_block.LOAD_BLOCK)
    stage = simulation[1]
    @test stage[:x_storage]["Lake_Taupo"].in == 750.275
    @test sort(collect(axes(stage[:u_thermal], 1))) == sort(G)
    @test sort(collect(axes(stage[:price], 1))) == sort(P)
    @test sort(collect(axes(stage[:marginal_water_value], 1))) == sort(R)
    for stage in simulation
        # Prices and water values are in $/MWh, and non-negative.
        @test all(p -> isfinite(p) && p >= -1e-6, stage[:price])
        @test all(v -> isfinite(v) && v >= -1e-6, stage[:marginal_water_value])
        # Thermal generation is within capacity.
        for g in G, b in B
            station = only(filter(r -> r.GENERATOR == g, data.thermal_stations))
            @test -1e-6 <= stage[:u_thermal][g, b] <= station.CAPACITY + 1e-6
        end
        # Storage is carried between weeks.
        @test all(r -> 0 <= stage[:x_storage][r].out, R)
    end
    for t in 2:52
        for r in R
            @test simulation[t][:x_storage][r].in ≈
                  simulation[t-1][:x_storage][r].out
        end
    end
    return
end

function test_create_graph()
    @test_throws(
        AssertionError,
        JADE._create_graph(; start_week = 0, discount_factor = 0.9),
    )
    @test_throws(
        AssertionError,
        JADE._create_graph(; start_week = 53, discount_factor = 0.9),
    )
    @test_throws(
        AssertionError,
        JADE._create_graph(; start_week = 1, discount_factor = 1.0),
    )
    @test_throws(
        AssertionError,
        JADE._create_graph(; start_week = 1, discount_factor = 0.0),
    )
    graph = JADE._create_graph(; start_week = 1, discount_factor = 0.9)
    @test graph.root_node == 0
    @test length(graph.nodes) == 53
    @test graph.nodes[0] == [(1, 1.0)]
    @test graph.nodes[1] == [(2, 1.0)]
    @test graph.nodes[22] == [(23, 1.0)]
    @test graph.nodes[23] == [(24, 1.0)]
    @test graph.nodes[52] == [(1, 0.9)]
    graph = JADE._create_graph(; start_week = 23, discount_factor = 0.95)
    @test graph.root_node == 0
    @test length(graph.nodes) == 53
    @test graph.nodes[0] == [(23, 1.0)]
    @test graph.nodes[1] == [(2, 1.0)]
    @test graph.nodes[22] == [(23, 0.95)]
    @test graph.nodes[23] == [(24, 1.0)]
    @test graph.nodes[52] == [(1, 1.0)]
    return
end

end  # TestJADE

TestJADE.runtests()
