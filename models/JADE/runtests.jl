# Copyright (c) 2016-26: Electric Power Optimization Centre, University of
# Auckland, and Oscar Dowson and contributors.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

module TestJADE

using Test

import DataFrames
import DBInterface
import HiGHS
import JuMP
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
    @test size(data.inflows) == (26 * 48 * 52, 4)
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
    @test size(df, 1) == 26 * 48 * 52
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
# Write one DataFrame per table to a temporary database and read it back.
function _data_from_tables(; kwargs...)
    filename = joinpath(mktempdir(), "data.db")
    db = SQLite.DB(filename)
    JADE._create_tables(db)
    for table in JADE.TABLES
        SQLite.load!(kwargs[Symbol(table)], db, table)
    end
    DBInterface.close!(db)
    return JADE.Data(filename)
end

# The tables of `data`, with some replaced by the keyword arguments.
function _data_with_tables(data::JADE.Data; kwargs...)
    tables = Dict(
        Symbol(t) => get(kwargs, Symbol(t), getfield(data, Symbol(t))) for
        t in JADE.TABLES
    )
    return _data_from_tables(; tables...)
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
    @test length(model[1].noise_terms) == 48
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
    for inflow_years in ([1995], 1980:1984, [2001, 1973, 2017, 1973])
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
    # 1970 to 2016: the last year cannot start a scenario
    @test length(scenarios) == 47
    for scenario in scenarios
        @test first.(scenario) == [39:52; 1:38]
    end
    # The first scenario is week 39 of 1970 to week 38 of 1971
    @test scenarios[1][1][2] == _historical_inflows(data, 1970, 39)
    @test scenarios[1][14][2] == _historical_inflows(data, 1970, 52)
    @test scenarios[1][15][2] == _historical_inflows(data, 1971, 1)
    @test scenarios[1][end][2] == _historical_inflows(data, 1971, 38)
    # The last scenario is week 39 of 2016 to week 38 of 2017
    @test scenarios[end][1][2] == _historical_inflows(data, 2016, 39)
    @test scenarios[end][end][2] == _historical_inflows(data, 2017, 38)
    # If the model starts in week 1, every year is a scenario
    config_1 = JADE.Config(; optimizer = HiGHS.Optimizer)
    sampling_scheme = JADE.HistoricalYears(data, config_1)
    scenarios = [noise.term for noise in sampling_scheme.scenarios]
    @test length(scenarios) == 48
    for scenario in scenarios
        @test first.(scenario) == 1:52
    end
    @test scenarios[1][1][2] == _historical_inflows(data, 1970, 1)
    @test scenarios[end][end][2] == _historical_inflows(data, 2017, 52)
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
    # Week 39 of 1980 to week 38 of 2017
    @test length(scenario) == 52 * (2017 - 1980)
    @test scenario[1] == (39, _historical_inflows(data, 1980, 39))
    @test scenario[14] == (52, _historical_inflows(data, 1980, 52))
    @test scenario[15] == (1, _historical_inflows(data, 1981, 1))
    @test scenario[end] == (38, _historical_inflows(data, 2017, 38))
    # The default start year is the first historical year
    sampling_scheme = JADE.HistoricalSequence(data, config)
    scenario_default = only(sampling_scheme.scenarios).term
    @test length(scenario_default) == 52 * (2017 - 1970)
    @test scenario_default[1] == (39, _historical_inflows(data, 1970, 39))
    @test scenario_default[end] == (38, _historical_inflows(data, 2017, 38))
    # If the model starts in week 1, the scenario ends in week 52
    config_1 = JADE.Config(; optimizer = HiGHS.Optimizer)
    sampling_scheme = JADE.HistoricalSequence(data, config_1; start_year = 2016)
    scenario_1 = only(sampling_scheme.scenarios).term
    @test length(scenario_1) == 52 * 2
    @test scenario_1[1] == (1, _historical_inflows(data, 2016, 1))
    @test scenario_1[end] == (52, _historical_inflows(data, 2017, 52))
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
        start_year = 2016,
        end_year = 2016,
    )
    scenario_end_1 = only(sampling_scheme.scenarios).term
    @test length(scenario_end_1) == 52
    @test scenario_end_1[1] == (1, _historical_inflows(data, 2016, 1))
    @test scenario_end_1[end] == (52, _historical_inflows(data, 2016, 52))
    # Invalid start and end years
    @test_throws(
        AssertionError,
        JADE.HistoricalSequence(data, config; start_year = 1900),
    )
    @test_throws(
        AssertionError,
        JADE.HistoricalSequence(data, config; end_year = 2020),
    )
    # Week 39 of 2017 would need inflows from 2018
    @test_throws(
        AssertionError,
        JADE.HistoricalSequence(data, config; start_year = 2017),
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
    sampling_scheme = JADE.HistoricalSequence(data, config; start_year = 2017)
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

# ## The toy system
#
# The default data is too large to check the subproblem by hand, so these tests
# use a toy system with closed-form solutions:
#
#  * Two nodes, `A` and `B`, joined by a transmission line from `A` to `B`
#    with `tranches` of `(capacity, loss_fraction)`.
#  * One thermal station `T` at `A`. Its cost is
#    OMCOST + HEAT_RATE * (fuel price + CO2 content * CO2 price)
#    = 5 + 10 * (2 + 0.05 * 20) = 35 $/MWh.
#  * One fixed station `F` at `B` with `fixed_output` MW in every block.
#  * All of the demand is at `B`: `demand` MW in each block. Lost load costs
#    1000 $/MWh at both nodes.
#  * A hydro network with a reservoir `R`, a junction `J`, and a hydro station
#    `H` at `A` from `R` to `J` with a capacity of 50 MW and a specific power
#    of 1 MW/cumec. Water spills from `R` to `J`, and `J` drains to the `SEA`
#    with a minimum flow of `min_flow`. `R` starts at `initial_level` Mm³ with
#    a maximum of 100 Mm³, and has an `inflow` of cumecs in every week.
#  * Two load blocks, `PEAK` with 40 hours and `OFF` with 128 hours.
#  * An `outage` of `T` in week 2.

const TOY_BLOCKS = ["PEAK", "OFF"]
const TOY_HOURS = Dict("PEAK" => 40.0, "OFF" => 128.0)
const TOY_THERMAL_COST = 35.0
const TOY_LOST_LOAD_COST = 1_000.0

function _toy_data(;
    demand::Dict{String,Float64} = Dict("PEAK" => 100.0, "OFF" => 60.0),
    tranches::Vector{Tuple{Float64,Float64}} = [(1_000.0, 0.1)],
    thermal_capacity::Float64 = 1_000.0,
    fixed_output::Float64 = 0.0,
    outage::Float64 = 0.0,
    initial_level::Float64 = 0.0,
    inflow::Float64 = 0.0,
    min_flow::Float64 = 0.0,
)
    DataFrame = DataFrames.DataFrame
    weeks, B = 1:52, TOY_BLOCKS
    n_outages = outage > 0 ? 1 : 0
    return _data_from_tables(;
        hours_per_block = DataFrame([
            (WEEK = w, LOAD_BLOCK = b, HOURS = TOY_HOURS[b]) for w in weeks for
            b in B
        ],),
        demand = DataFrame([
            (
                NODE = n,
                WEEK = w,
                LOAD_BLOCK = b,
                DEMAND = n == "B" ? demand[b] : 0.0,
            ) for n in ["A", "B"] for w in weeks for b in B
        ],),
        lost_load = DataFrame(;
            NODE = ["A", "B"],
            COST = [TOY_LOST_LOAD_COST, TOY_LOST_LOAD_COST],
            PROPORTION = [1.0, 1.0],
        ),
        fixed_stations = DataFrame(; STATION = ["F"], NODE = ["B"]),
        fixed_station_output = DataFrame([
            (STATION = "F", LOAD_BLOCK = b, OUTPUT = fixed_output) for b in B
        ],),
        thermal_stations = DataFrame(;
            GENERATOR = ["T"],
            NODE = ["A"],
            FUEL = ["gas"],
            HEAT_RATE = [10.0],
            CAPACITY = [thermal_capacity],
            OMCOST = [5.0],
        ),
        fuels = DataFrame(; FUEL = ["gas"], CO2_CONTENT = [0.05]),
        fuel_prices = DataFrame([
            (WEEK = w, FUEL = "gas", PRICE = 2.0) for w in weeks
        ],),
        co2_prices = DataFrame(;
            WEEK = collect(weeks),
            CO2_PRICE = fill(20.0, 52),
        ),
        hydro_stations = DataFrame(;
            GENERATOR = ["H"],
            HEAD_WATER_FROM = ["R"],
            TAIL_WATER_TO = ["J"],
            NODE = ["A"],
            CAPACITY = [50.0],
            SPECIFIC_POWER = [1.0],
        ),
        hydro_arcs = DataFrame(;
            ORIG = ["R", "J"],
            DEST = ["J", "SEA"],
            MIN_FLOW = [0.0, min_flow],
            MAX_FLOW = [missing, missing],
        ),
        reservoirs = DataFrame(;
            RESERVOIR = ["R"],
            INITIAL_LEVEL = [initial_level],
            MAX_LEVEL = [100.0],
        ),
        inflows = DataFrame([
            (CATCHMENT = "R", YEAR = 2000, WEEK = w, INFLOW = inflow) for
            w in weeks
        ],),
        station_outages = DataFrame(;
            GENERATOR = fill("T", n_outages),
            WEEK = fill(2, n_outages),
            OUTAGE = fill(outage, n_outages),
        ),
        transmission = DataFrame([
            (
                FROM_NODE = "A",
                TO_NODE = "B",
                TRANCHE = i,
                CAPACITY = c,
                LOSS_FRACTION = l,
            ) for (i, (c, l)) in enumerate(tranches)
        ],),
    )
end

# Simulate the first two weeks of the toy system with an untrained policy. The
# cost-to-go is the lower bound of zero, so each week is solved on its own and
# the solution is the cheapest way to meet this week's demand.
function _simulate_toy(data::JADE.Data)
    config = JADE.Config(; optimizer = HiGHS.Optimizer)
    model = SDDP.PolicyGraph(data, config)
    simulations = SDDP.simulate(
        model,
        1,
        [
            :u_thermal,
            :u_transmission,
            :u_lost_load,
            :x_storage,
            :u_turbine,
            :u_flow,
            :u_flow_below,
        ];
        sampling_scheme = SDDP.InSampleMonteCarlo(;
            max_depth = 2,
            terminate_on_dummy_leaf = false,
        ),
        custom_recorders = Dict{Symbol,Function}(
            :price => JADE._price,
            :water_value => JADE._marginal_water_value,
        ),
    )
    return model, only(simulations)
end

_toy_line(i::Int = 1) = ("A", "B", i)

function test_toy_c_demand()
    demand = Dict("PEAK" => 100.0, "OFF" => 60.0)
    data = _toy_data(; demand, fixed_output = 10.0)
    config = JADE.Config(; optimizer = HiGHS.Optimizer)
    model = SDDP.PolicyGraph(data, config)
    sp = model[1].subproblem
    line = _toy_line()
    for b in TOY_BLOCKS
        # Half of the losses are attributed to each end of the line: the
        # sending node pays (1 + 0.5 * loss) and the receiving node gets
        # (1 - 0.5 * loss).
        c_A, c_B = sp[:c_demand]["A", b], sp[:c_demand]["B", b]
        u_line = sp[:u_transmission][line, b]
        @test JuMP.normalized_coefficient(c_A, u_line) ≈ -(1 + 0.5 * 0.1)
        @test JuMP.normalized_coefficient(c_B, u_line) ≈ 1 - 0.5 * 0.1
        # Thermal generation and lost load count only at their own node.
        u_thermal = sp[:u_thermal]["T", b]
        @test JuMP.normalized_coefficient(c_A, u_thermal) == 1.0
        @test JuMP.normalized_coefficient(c_B, u_thermal) == 0.0
        u_lost_load = sp[:u_lost_load][("B", TOY_LOST_LOAD_COST), b]
        @test JuMP.normalized_coefficient(c_A, u_lost_load) == 0.0
        @test JuMP.normalized_coefficient(c_B, u_lost_load) == 1.0
        # Hydro generation is specific power times the turbine flow.
        u_turbine = sp[:u_turbine]["H", b]
        @test JuMP.normalized_coefficient(c_A, u_turbine) == 1.0
        @test JuMP.normalized_coefficient(c_B, u_turbine) == 0.0
        # Fixed generation is a constant, so it moves to the right-hand side.
        @test JuMP.normalized_rhs(c_A) == 0.0
        @test JuMP.normalized_rhs(c_B) == demand[b] - 10.0
        @test JuMP.constraint_object(c_B).set isa JuMP.MOI.GreaterThan
    end
    return
end

function test_toy_transmission_losses()
    loss = 0.1
    data = _toy_data(; tranches = [(1_000.0, loss)])
    _, simulation = _simulate_toy(data)
    week = simulation[1]
    for (b, demand) in Dict("PEAK" => 100.0, "OFF" => 60.0)
        # The flow is the demand plus the receiving half of the losses, and
        # the thermal generation is the flow plus the sending half.
        flow = demand / (1 - 0.5 * loss)
        @test week[:u_transmission][_toy_line(), b] ≈ flow atol = 1e-6
        @test week[:u_thermal]["T", b] ≈ (1 + 0.5 * loss) * flow atol = 1e-6
        # The total losses are the loss fraction of the flow.
        @test week[:u_thermal]["T", b] - demand ≈ loss * flow atol = 1e-6
        @test week[:u_lost_load][("B", TOY_LOST_LOAD_COST), b] ≈ 0 atol = 1e-6
        # The price at `A` is the thermal cost. The price at `B` is marked up
        # by the losses.
        @test week[:price]["A", b] ≈ TOY_THERMAL_COST atol = 1e-6
        @test week[:price]["B", b] ≈
              TOY_THERMAL_COST * (1 + 0.5 * loss) / (1 - 0.5 * loss) atol = 1e-6
    end
    # The stage objective is the scaled cost of the thermal generation.
    cost = sum(
        TOY_HOURS[b] * TOY_THERMAL_COST * week[:u_thermal]["T", b] for
        b in TOY_BLOCKS
    )
    @test week[:stage_objective] ≈ 1e-6 * cost atol = 1e-9
    return
end

function test_toy_transmission_tranches()
    tranches = [(50.0, 0.1), (1_000.0, 0.3)]
    demand = Dict("PEAK" => 100.0, "OFF" => 40.0)
    data = _toy_data(; tranches, demand)
    _, simulation = _simulate_toy(data)
    week = simulation[1]
    # In `PEAK`, the first tranche is full and the second carries the rest.
    u_1, u_2 = (week[:u_transmission][_toy_line(i), "PEAK"] for i in 1:2)
    @test u_1 ≈ 50.0 atol = 1e-6
    @test u_2 ≈ (100.0 - (1 - 0.05) * 50.0) / (1 - 0.15) atol = 1e-6
    @test week[:u_thermal]["T", "PEAK"] ≈ 1.05 * u_1 + 1.15 * u_2 atol = 1e-6
    @test week[:price]["B", "PEAK"] ≈ TOY_THERMAL_COST * 1.15 / 0.85 atol = 1e-6
    # In `OFF`, the demand fits in the first tranche.
    u_1, u_2 = (week[:u_transmission][_toy_line(i), "OFF"] for i in 1:2)
    @test u_1 ≈ 40.0 / 0.95 atol = 1e-6
    @test u_2 ≈ 0.0 atol = 1e-6
    @test week[:price]["B", "OFF"] ≈ TOY_THERMAL_COST * 1.05 / 0.95 atol = 1e-6
    return
end

function test_toy_lost_load()
    data = _toy_data(; thermal_capacity = 60.0)
    _, simulation = _simulate_toy(data)
    week = simulation[1]
    # Every MW sent from `A` delivers (1 - 0.05) / (1 + 0.05) MW at `B`.
    delivered = 60.0 * 0.95 / 1.05
    for (b, demand) in Dict("PEAK" => 100.0, "OFF" => 60.0)
        @test week[:u_thermal]["T", b] ≈ 60.0 atol = 1e-6
        lost_load = week[:u_lost_load][("B", TOY_LOST_LOAD_COST), b]
        @test lost_load ≈ demand - delivered atol = 1e-6
        # The price at `B` is the cost of lost load, and the price at `A` is
        # the value of the load that one more MW at `A` would serve.
        @test week[:price]["B", b] ≈ TOY_LOST_LOAD_COST atol = 1e-6
        @test week[:price]["A", b] ≈ TOY_LOST_LOAD_COST * 0.95 / 1.05 atol =
            1e-6
    end
    cost = sum(
        TOY_HOURS[b] * (
            TOY_THERMAL_COST * week[:u_thermal]["T", b] +
            TOY_LOST_LOAD_COST *
            week[:u_lost_load][("B", TOY_LOST_LOAD_COST), b]
        ) for b in TOY_BLOCKS
    )
    @test week[:stage_objective] ≈ 1e-6 * cost atol = 1e-9
    return
end

function test_toy_outage()
    data = _toy_data(; outage = 950.0)
    model, simulation = _simulate_toy(data)
    # The outage is in week 2 only.
    @test JuMP.upper_bound(model[1].subproblem[:u_thermal]["T", "PEAK"]) ==
          1_000.0
    @test JuMP.upper_bound(model[2].subproblem[:u_thermal]["T", "PEAK"]) == 50.0
    @test simulation[1][:u_thermal]["T", "PEAK"] ≈ 100.0 / 0.95 * 1.05 atol =
        1e-6
    @test simulation[2][:u_thermal]["T", "PEAK"] ≈ 50.0 atol = 1e-6
    @test simulation[2][:u_lost_load][("B", TOY_LOST_LOAD_COST), "PEAK"] ≈
          100.0 - 50.0 * 0.95 / 1.05 atol = 1e-6
    return
end

function test_toy_free_disposal()
    # The fixed station exceeds the demand, so the surplus is disposed of
    # without cost: `c_demand` is `>=`, not `==`.
    data = _toy_data(; fixed_output = 200.0)
    _, simulation = _simulate_toy(data)
    week = simulation[1]
    for b in TOY_BLOCKS
        @test week[:u_thermal]["T", b] ≈ 0.0 atol = 1e-6
        @test week[:u_transmission][_toy_line(), b] ≈ 0.0 atol = 1e-6
        @test week[:u_lost_load][("B", TOY_LOST_LOAD_COST), b] ≈ 0.0 atol = 1e-6
        @test week[:price]["B", b] ≈ 0.0 atol = 1e-6
    end
    @test week[:stage_objective] ≈ 0.0 atol = 1e-9
    return
end

function test_toy_hydro_storage()
    # With 10 Mm³ of water and no inflow, the untrained policy uses all of
    # the water this week to displace thermal generation.
    data = _toy_data(; initial_level = 10.0)
    _, simulation = _simulate_toy(data)
    week = simulation[1]
    @test week[:x_storage]["R"].in == 10.0
    @test week[:x_storage]["R"].out ≈ 0.0 atol = 1e-6
    hydro_energy =
        sum(TOY_HOURS[b] * week[:u_turbine]["H", b] for b in TOY_BLOCKS)
    @test hydro_energy ≈ 10.0 / JADE.MM3_PER_CUMEC_HOUR atol = 1e-6
    # Every MWh of hydro displaces a MWh of thermal, so the marginal value of
    # water is the thermal cost.
    @test week[:water_value]["R"] ≈ TOY_THERMAL_COST atol = 1e-6
    for b in TOY_BLOCKS
        @test week[:price]["A", b] ≈ TOY_THERMAL_COST atol = 1e-6
    end
    # With 100 Mm³, the capacity of the station binds instead, and the water
    # has no value because the untrained policy has no future.
    data = _toy_data(; initial_level = 100.0)
    _, simulation = _simulate_toy(data)
    week = simulation[1]
    for b in TOY_BLOCKS
        @test week[:u_turbine]["H", b] ≈ 50.0 atol = 1e-6
    end
    water_used = JADE.MM3_PER_CUMEC_HOUR * 50.0 * sum(values(TOY_HOURS))
    @test week[:x_storage]["R"].out ≈ 100.0 - water_used atol = 1e-6
    @test week[:water_value]["R"] ≈ 0.0 atol = 1e-6
    # The state is passed to the next week.
    @test simulation[2][:x_storage]["R"].in ≈ 100.0 - water_used atol = 1e-6
    return
end

function test_toy_hydro_spill()
    # The reservoir is full and the inflow exceeds the capacity of the
    # station, so the excess must spill through the junction to the sea.
    # Spilled water is free, so the solution is degenerate: the policy may
    # also drain the reservoir and spill the water. We check the balances
    # rather than a particular solution.
    data = _toy_data(; initial_level = 100.0, inflow = 100.0)
    _, simulation = _simulate_toy(data)
    week = simulation[1]
    storage_in, storage_out =
        week[:x_storage]["R"].in, week[:x_storage]["R"].out
    @test storage_in == 100.0
    @test -1e-6 <= storage_out <= 100.0 + 1e-6
    net_flow = 0.0
    for b in TOY_BLOCKS
        turbine, spill = week[:u_turbine]["H", b], week[:u_flow][("R", "J"), b]
        @test turbine ≈ 50.0 atol = 1e-6
        @test spill >= -1e-6
        # The junction cannot store water.
        @test week[:u_flow][("J", "SEA"), b] ≈ turbine + spill atol = 1e-6
        net_flow +=
            JADE.MM3_PER_CUMEC_HOUR * TOY_HOURS[b] * (100.0 - turbine - spill)
    end
    # The reservoir balance converts cumecs to Mm³ over the hours of each block.
    @test storage_out - storage_in ≈ net_flow atol = 1e-6
    @test week[:water_value]["R"] ≈ 0.0 atol = 1e-6
    return
end

function test_toy_min_flow_penalty()
    # Only 100 cumecs of water reach the junction over the week, so a minimum
    # flow of 150 cumecs is violated by 50 cumecs on average. The reservoir is
    # balanced over the week, not within it, so the policy may move water
    # between the blocks; only the total violation in cumec-hours is fixed.
    # The penalty is converted from $/MWh with the largest specific power in
    # the system, 1 MW/cumec.
    data = _toy_data(; inflow = 100.0, min_flow = 150.0)
    _, simulation = _simulate_toy(data)
    week = simulation[1]
    hours = sum(values(TOY_HOURS))
    for b in TOY_BLOCKS
        flow = week[:u_flow][("J", "SEA"), b]
        @test flow + week[:u_flow_below][("J", "SEA"), b] ≈ 150.0 atol = 1e-6
    end
    violation = sum(
        TOY_HOURS[b] * week[:u_flow_below][("J", "SEA"), b] for b in TOY_BLOCKS
    )
    @test violation ≈ 50.0 * hours atol = 1e-6
    @test week[:x_storage]["R"].out ≈ 0.0 atol = 1e-6
    thermal_cost = sum(
        TOY_HOURS[b] * TOY_THERMAL_COST * week[:u_thermal]["T", b] for
        b in TOY_BLOCKS
    )
    penalty = 500.0 * 1.0 * violation
    @test week[:stage_objective] ≈ 1e-6 * (thermal_cost + penalty) atol = 1e-9
    return
end

end  # TestJADE

TestJADE.runtests()
