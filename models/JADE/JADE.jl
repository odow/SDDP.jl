# Copyright (c) 2016-26: Electric Power Optimization Centre, University of #src
# Auckland, and Oscar Dowson and contributors.                             #src
#                                                                          #src
# This Source Code Form is subject to the terms of the Mozilla Public      #src
# License, v2.0. If a copy of the MPL was not distributed with this file,  #src
# You can obtain one at http://mozilla.org/MPL/2.0/.                       #src

# # JADE
#
# This file is a simplified copy of
# [JADE.jl](https://github.com/EPOC-NZ/JADE.jl), the New Zealand Electricity
# Authority's model for the medium-term hydro-thermal scheduling problem. JADE
# stands for Just Another DOASA Environment.
#
# It builds the model as an infinite-horizon
# [SDDP.jl](https://github.com/odow/SDDP.jl) policy graph with 52 weekly nodes,
# where the uncertainty is the weekly inflow into the hydro catchments. This
# file only builds the model. Training, simulating, and analysing the results
# are done with SDDP.jl.
#
# The source code is a single file, `models/JADE/JADE.jl`, in the
# [SDDP.jl repository](https://github.com/odow/SDDP.jl/tree/master/models/JADE)
# on GitHub. Download it, along with the `data` directory beside it, to use or
# modify it. You can also find a local copy in the SDDP.jl package:
# ```julia
# using SDDP
# JADE_DIR = joinpath(pkgdir(SDDP), "models", "JADE")
# ```
#
# ## Purpose
#
# This file has two purposes:
#
#  1. It is a working, minimum viable implementation of JADE.
#  2. It is a guide to the structure and layout of a good SDDP model. If you
#     are building your own SDDP.jl model, read this file for an example of how
#     to lay out the data, the policy graph, the subproblem, and the sampling
#     schemes.
#
# The file is written as a literate program. You should be able to start at the
# top of this file and read to the bottom, and understand everything that is
# contained within it.
#
# A single file is also easy to modify. We do not provide options for every
# way that the model could be extended. Instead, if you want to extend the
# model, for example, by adding a new type of generator, a new constraint, or a
# different graph, copy this file and change it. The file is short enough that
# you can give it to an AI agent and ask it to make the change for you.
#
# Non-goals:
#
#  * This file is not a flexible implementation of an arbitrary hydro-thermal
#    scheduling model. It is opinionated. Where there is a choice, we make it
#    rather than adding an option.
#  * The input data is the minimum required data to implement the model. The
#    purpose is not to provide a flexible post-hoc analysis tool for different
#    data sets; we leave the choice of that to the user.
#
# ## Usage
#
# The input to the model is a single SQLite database. The intended usage is:
#
# ```julia
# using SDDP
# import HiGHS
# JADE_DIR = joinpath(pkgdir(SDDP), "models", "JADE")
# include(joinpath(JADE_DIR, "JADE.jl"))
# data = JADE.Data(;
#     directory = joinpath(JADE_DIR, "data"),
#     filename = joinpath(JADE_DIR, "data.db"),
# )
# config = JADE.Config(; optimizer = HiGHS.Optimizer)
# model = SDDP.PolicyGraph(data, config)
# JADE.train(model; iteration_limit = 100)
# simulations = JADE.simulate(model, 100)
# ```
#
# `JADE.train` wraps `SDDP.train` to choose a sampling scheme that suits the
# cyclic policy graph (see the section "Training the policy" below).
#
# `JADE.simulate` wraps `SDDP.simulate` to record the reservoir storage, the
# thermal generation, the price of energy at each node, and the marginal value
# of water in each reservoir.
#
# `JADE.Config` holds the options that a user can change, such as the week of
# the year in which the model starts, or the historical years used as the
# in-sample inflows. See the docstring of `Config` below.
#
# To evaluate a policy on the historical inflows, rather than on inflows sampled
# independently in each week, use `JADE.HistoricalYears`, which has one
# scenario for each historical year, or `JADE.HistoricalSequence`, which has a
# single scenario that runs through the historical record:
#
# ```julia
# sampling_scheme = SDDP.JADE.HistoricalYears(data, config)
# n = length(sampling_scheme.scenarios)
# simulations = JADE.simulate(model, n; sampling_scheme)
# ```
#
# If your data is in a directory of CSV files, pass the directory as well, and
# `JADE.Data` converts it to a database at `filename` before reading it (see the
# last section of this file). Pass `force = true` to overwrite an existing file:
#
# ```julia
# data = JADE.Data(; filename = "data.db", directory = "data", force = true)
# ```
#
# ## The module
#
# Everything is defined in the `JADE` module. JuMP and SDDP.jl build the
# model, SQLite.jl and DataFrames.jl read the input data, and CSV.jl reads the
# CSV files that can be converted to the input data.

module JADE

import CSV
import DataFrames
import DBInterface
import JuMP
import SDDP
import SQLite

# ## The input data
#
# This section defines the input data, and how it is read and validated.
#
# The input data is a SQLite database with the following tables:
#
# | Table                  | Contents                                         |
# | :--------------------- | :----------------------------------------------- |
# | `hours_per_block`      | Hours in each load block of each week            |
# | `demand`               | Demand at each node in each week and load block  |
# | `lost_load`            | Sheddable demand tranches and their cost         |
# | `fixed_stations`       | Node of each fixed-output (e.g., wind) station   |
# | `fixed_station_output` | Output of each fixed station in each load block  |
# | `thermal_stations`     | Node, fuel, heat rate, capacity, and O&M cost    |
# | `fuels`                | CO₂ content of each fuel                         |
# | `fuel_prices`          | Price of each fuel in each week                  |
# | `co2_prices`           | Price of CO₂ in each week                        |
# | `hydro_stations`       | Turbines between two hydro nodes                 |
# | `hydro_arcs`           | Rivers, canals, and spillways; min and max flows |
# | `reservoirs`           | Initial and maximum storage                      |
# | `inflows`              | Historical inflows by catchment, year, and week  |
# | `station_outages`      | Reductions in station capacity in each week      |
# | `transmission`         | Directed lines, with piecewise-linear losses     |
#
# `data.db` is a small three-node (`NI`, `HAY`, `SI`) model of New Zealand with
# 48 years (1970–2017) of historical inflows.
#
# ### The schema
#
# This schema is the single source of truth for the input data: every column,
# its units, and the simple constraints on its values are defined here, and
# SQLite enforces them when the data is written. Tables are `STRICT`, so SQLite
# also enforces the type of each column. The tables are in third normal form:
# every table has a primary key, and every other column depends on the whole
# key and nothing else. Tables are `WITHOUT ROWID`, which makes the file
# smaller.
#
# Every table uses the same conventions:
#
#  * Time is indexed by `WEEK`, where week `1` is the first week of January, and
#    there are exactly 52 weeks in a year. There are no years (except in the
#    historical inflow record); every weekly time series repeats each year.
#  * The load blocks, like `PEAK` or `OFFPEAK`, are arbitrary names, but they
#    must be the same in every table.
#  * `NULL` is used only to mean "there is no limit."
#  * The node `SEA` is special. It is the sink of the hydro network.

const SCHEMA = """
-- The number of hours in each load block of each week.
CREATE TABLE hours_per_block (
    WEEK INTEGER NOT NULL CHECK (WEEK BETWEEN 1 AND 52),
    LOAD_BLOCK TEXT NOT NULL,
    HOURS REAL NOT NULL CHECK (HOURS >= 0),        -- hours
    PRIMARY KEY (WEEK, LOAD_BLOCK)
) STRICT, WITHOUT ROWID;

-- The average demand at each node of the power network in each load block.
CREATE TABLE demand (
    NODE TEXT NOT NULL,
    WEEK INTEGER NOT NULL CHECK (WEEK BETWEEN 1 AND 52),
    LOAD_BLOCK TEXT NOT NULL,
    DEMAND REAL NOT NULL CHECK (DEMAND >= 0),      -- MW
    PRIMARY KEY (NODE, WEEK, LOAD_BLOCK)
) STRICT, WITHOUT ROWID;

-- Tranches of demand at each node that can be shed at a cost. The size of a
-- tranche is PROPORTION of the demand in each load block. The tranches of a
-- node are identified by their COST.
CREATE TABLE lost_load (
    NODE TEXT NOT NULL,
    COST REAL NOT NULL CHECK (COST >= 0),          -- \$/MWh
    PROPORTION REAL NOT NULL CHECK (PROPORTION BETWEEN 0 AND 1),
    PRIMARY KEY (NODE, COST)
) STRICT, WITHOUT ROWID;

-- Generators, typically wind, whose output is fixed. Their output is
-- subtracted from demand.
CREATE TABLE fixed_stations (
    STATION TEXT NOT NULL PRIMARY KEY,
    NODE TEXT NOT NULL
) STRICT, WITHOUT ROWID;

-- The output of each fixed station in each load block. It is the same in every
-- week.
CREATE TABLE fixed_station_output (
    STATION TEXT NOT NULL,
    LOAD_BLOCK TEXT NOT NULL,
    OUTPUT REAL NOT NULL CHECK (OUTPUT >= 0),      -- MW
    PRIMARY KEY (STATION, LOAD_BLOCK)
) STRICT, WITHOUT ROWID;

-- Thermal generators.
CREATE TABLE thermal_stations (
    GENERATOR TEXT NOT NULL PRIMARY KEY,
    NODE TEXT NOT NULL,
    FUEL TEXT NOT NULL,
    HEAT_RATE REAL NOT NULL CHECK (HEAT_RATE >= 0),  -- GJ/MWh
    CAPACITY REAL NOT NULL CHECK (CAPACITY >= 0),    -- MW
    OMCOST REAL NOT NULL CHECK (OMCOST >= 0)         -- \$/MWh
) STRICT, WITHOUT ROWID;

-- The fuels used by thermal stations, and the carbon emitted when they are
-- burnt.
CREATE TABLE fuels (
    FUEL TEXT NOT NULL PRIMARY KEY,
    CO2_CONTENT REAL NOT NULL CHECK (CO2_CONTENT >= 0)   -- tonnes CO2/GJ
) STRICT, WITHOUT ROWID;

-- The price of each fuel in each week. The total cost of a fuel in a week is
-- PRICE + CO2_CONTENT * CO2_PRICE in \$/GJ.
CREATE TABLE fuel_prices (
    WEEK INTEGER NOT NULL CHECK (WEEK BETWEEN 1 AND 52),
    FUEL TEXT NOT NULL,
    PRICE REAL NOT NULL CHECK (PRICE >= 0),              -- \$/GJ
    PRIMARY KEY (WEEK, FUEL)
) STRICT, WITHOUT ROWID;

-- The price of carbon emissions in each week.
CREATE TABLE co2_prices (
    WEEK INTEGER NOT NULL PRIMARY KEY CHECK (WEEK BETWEEN 1 AND 52),
    CO2_PRICE REAL NOT NULL CHECK (CO2_PRICE >= 0)       -- \$/tonne CO2
) STRICT, WITHOUT ROWID;

-- Hydro generators. Each row is a turbine that takes water from
-- HEAD_WATER_FROM to TAIL_WATER_TO and produces power at NODE. A spillway is
-- not part of a station: it is a row in hydro_arcs.
CREATE TABLE hydro_stations (
    GENERATOR TEXT NOT NULL PRIMARY KEY,
    HEAD_WATER_FROM TEXT NOT NULL,
    TAIL_WATER_TO TEXT NOT NULL,
    NODE TEXT NOT NULL,
    CAPACITY REAL NOT NULL CHECK (CAPACITY >= 0),              -- MW
    SPECIFIC_POWER REAL NOT NULL CHECK (SPECIFIC_POWER >= 0)   -- MW/cumec
) STRICT, WITHOUT ROWID;

-- Channels through which water can flow without generating power, including
-- rivers, canals, and spillways. A NULL MAX_FLOW means there is no limit.
CREATE TABLE hydro_arcs (
    ORIG TEXT NOT NULL,
    DEST TEXT NOT NULL,
    MIN_FLOW REAL NOT NULL CHECK (MIN_FLOW >= 0),   -- cumec
    MAX_FLOW REAL CHECK (MAX_FLOW >= MIN_FLOW),     -- cumec
    PRIMARY KEY (ORIG, DEST)
) STRICT, WITHOUT ROWID;

-- Nodes in the hydro network that can store water between weeks. The minimum
-- level is 0.
CREATE TABLE reservoirs (
    RESERVOIR TEXT NOT NULL PRIMARY KEY,
    INITIAL_LEVEL REAL NOT NULL CHECK (INITIAL_LEVEL >= 0),       -- Mm³
    MAX_LEVEL REAL NOT NULL CHECK (MAX_LEVEL >= INITIAL_LEVEL)    -- Mm³
) STRICT, WITHOUT ROWID;

-- The historical inflows into each node of the hydro network. Every CATCHMENT
-- must have an inflow in every week of every historical YEAR. In week t, the
-- uncertainty is sampled from week t of every YEAR.
CREATE TABLE inflows (
    CATCHMENT TEXT NOT NULL,
    YEAR INTEGER NOT NULL,
    WEEK INTEGER NOT NULL CHECK (WEEK BETWEEN 1 AND 52),
    INFLOW REAL NOT NULL,                           -- cumec
    PRIMARY KEY (CATCHMENT, YEAR, WEEK)
) STRICT, WITHOUT ROWID;

-- The reduction in the capacity of a hydro or thermal station in a week. Weeks
-- that are not listed have no outage.
CREATE TABLE station_outages (
    GENERATOR TEXT NOT NULL,
    WEEK INTEGER NOT NULL CHECK (WEEK BETWEEN 1 AND 52),
    OUTAGE REAL NOT NULL CHECK (OUTAGE >= 0),       -- MW
    PRIMARY KEY (GENERATOR, WEEK)
) STRICT, WITHOUT ROWID;

-- Directed transmission lines. A line in both directions requires two sets of
-- rows. The flow on a line is split into tranches, numbered from 1. Tranche k
-- carries up to CAPACITY MW and loses LOSS_FRACTION of the flow through it. The
-- capacity of the line is the sum of the capacities of its tranches.
CREATE TABLE transmission (
    FROM_NODE TEXT NOT NULL,
    TO_NODE TEXT NOT NULL,
    TRANCHE INTEGER NOT NULL CHECK (TRANCHE >= 1),
    CAPACITY REAL NOT NULL CHECK (CAPACITY > 0),                     -- MW
    LOSS_FRACTION REAL NOT NULL CHECK (LOSS_FRACTION BETWEEN 0 AND 1),
    PRIMARY KEY (FROM_NODE, TO_NODE, TRANCHE)
) STRICT, WITHOUT ROWID;
"""

# `TABLES` is the name of each table, in the order that it appears in the
# schema:

const TABLES = [String(m[1]) for m in eachmatch(r"CREATE TABLE (\w+)", SCHEMA)]

# `_create_tables` creates the tables in an empty database `db`. SQLite can
# execute only one statement at a time, so we split the schema on `;`. This
# means that the comments in the schema must not contain a `;`.

function _create_tables(db::SQLite.DB)
    for statement in split(SCHEMA, ";"; keepempty = false)
        if !isempty(strip(statement))
            DBInterface.execute(db, statement)
        end
    end
    return
end

# `_column_names` returns the names of the columns of `table` in `db`, in order:

function _column_names(db::SQLite.DB, table::String)
    info = DataFrames.DataFrame(
        DBInterface.execute(db, "PRAGMA table_info($table)"),
    )
    return Vector{String}(info.name)
end

# ### The `Data` type
#
# The tables are stored in a `Data` struct, with one field for each table in
# the schema. We define our own type, rather than using a `NamedTuple`, so that
# we can dispatch on it, for example, in `SDDP.PolicyGraph(data, config)`.
#
# The only way to create a `Data` is to read a database, so the struct has a
# single inner constructor, `Data(filename)`, and no constructor that takes the
# tables directly. The constructor reads every table of the database into a
# `DataFrames.DataFrame`. To give a helpful error if the database was not
# created from this schema, it first checks that the columns of each table match
# the columns of the same table in an empty in-memory database that was created
# from the schema. It then validates the data with `_validate_data`, which is
# defined in the next section.

"""
    Data(filename::String)
    Data(; filename::String, directory::String, force::Bool = false)

The input data: one `DataFrames.DataFrame` for each table in the schema. Access
a table as a field, for example, `data.reservoirs`.

`Data(filename)` reads the SQLite database `filename`.

`Data(; filename, directory)` first converts the CSV files in `directory` into
a new SQLite database `filename`, and then reads it. If `filename` already
exists, it is an error unless `force = true`, in which case it is overwritten.
"""
struct Data
    hours_per_block::DataFrames.DataFrame
    demand::DataFrames.DataFrame
    lost_load::DataFrames.DataFrame
    fixed_stations::DataFrames.DataFrame
    fixed_station_output::DataFrames.DataFrame
    thermal_stations::DataFrames.DataFrame
    fuels::DataFrames.DataFrame
    fuel_prices::DataFrames.DataFrame
    co2_prices::DataFrames.DataFrame
    hydro_stations::DataFrames.DataFrame
    hydro_arcs::DataFrames.DataFrame
    reservoirs::DataFrames.DataFrame
    inflows::DataFrames.DataFrame
    station_outages::DataFrames.DataFrame
    transmission::DataFrames.DataFrame

    function Data(filename::String)
        @assert isfile(filename) "Missing input file: $filename"
        expected = SQLite.DB()  # An empty in-memory database with the schema
        _create_tables(expected)
        db = SQLite.DB(filename)
        tables = try
            map(TABLES) do table
                @assert(
                    _column_names(db, table) == _column_names(expected, table),
                    "Table `$table` does not match the schema",
                )
                query = DBInterface.execute(db, "SELECT * FROM $table")
                return DataFrames.DataFrame(query)
            end
        finally
            DBInterface.close!(db)
        end
        data = new(tables...)  # The fields are in the order of `TABLES`
        _validate_data(data)
        return data
    end
end

# The fields must be in the same order as the tables in the schema:

@assert fieldnames(Data) == Tuple(Symbol.(TABLES))

# Printing a `Data` would otherwise print all fifteen tables.

function Base.show(io::IO, data::Data)
    print(io, "JADE.Data with $(length(TABLES)) tables")
    return
end

# ### Validating the data
#
# The schema checks each table on its own. Here, we check that the tables are
# consistent with each other. Catching these errors here gives a much better
# error message than a `KeyError` in the middle of building the model.

function _validate_data(data::Data)
    ## Every weekly time series must contain the weeks 1 to 52.
    weeks = 1:52
    @assert sort(unique(data.hours_per_block.WEEK)) == weeks
    @assert sort(unique(data.demand.WEEK)) == weeks
    @assert sort(unique(data.inflows.WEEK)) == weeks
    @assert sort(data.co2_prices.WEEK) == weeks
    ## Every table must use the same load blocks.
    blocks = sort(unique(data.hours_per_block.LOAD_BLOCK))
    @assert sort(unique(data.demand.LOAD_BLOCK)) == blocks
    @assert data.fixed_station_output.LOAD_BLOCK ⊆ blocks
    ## There must be one hours_per_block for each week and block, and one demand
    ## for each node, week, and block.
    @assert size(data.hours_per_block, 1) == length(weeks) * length(blocks)
    nodes = unique(data.demand.NODE)
    @assert size(data.demand, 1) ==
            length(nodes) * length(weeks) * length(blocks)
    ## Every fixed station must have an output in every block.
    stations = data.fixed_stations.STATION
    @assert data.fixed_station_output.STATION ⊆ stations
    @assert size(data.fixed_station_output, 1) ==
            length(stations) * length(blocks)
    ## Every node in the power system must have a demand.
    @assert data.fixed_stations.NODE ⊆ nodes
    @assert data.hydro_stations.NODE ⊆ nodes
    @assert data.thermal_stations.NODE ⊆ nodes
    @assert data.lost_load.NODE ⊆ nodes
    @assert data.transmission.FROM_NODE ⊆ nodes
    @assert data.transmission.TO_NODE ⊆ nodes
    ## Every node can shed all of its demand.
    for node in nodes
        tranches = filter(r -> r.NODE == node, data.lost_load)
        @assert sum(tranches.PROPORTION; init = 0.0) >= 1 "Lost load at $node"
    end
    ## Every fuel must have a price in every week.
    @assert data.thermal_stations.FUEL ⊆ data.fuels.FUEL
    @assert data.fuel_prices.FUEL ⊆ data.fuels.FUEL
    for fuel in data.fuels.FUEL
        fuel_weeks = filter(r -> r.FUEL == fuel, data.fuel_prices).WEEK
        @assert sort(fuel_weeks) == weeks "Fuel prices for $fuel"
    end
    ## Every outage must be of a known station, and no larger than its capacity.
    capacity = Dict(
        r.GENERATOR => r.CAPACITY for r in Iterators.flatten((
            eachrow(data.hydro_stations),
            eachrow(data.thermal_stations),
        ))
    )
    for r in eachrow(data.station_outages)
        @assert r.OUTAGE <= capacity[r.GENERATOR] "Outage of $(r.GENERATOR)"
    end
    ## The tranches of each transmission line must be numbered 1, 2, ..., and
    ## the loss fractions must be increasing so that the losses are convex.
    for line in _transmission_lines(data)
        losses = _transmission_losses(data, line)
        @assert issorted(last.(losses)) "Losses on $line must be increasing"
    end
    ## Every node in the hydro network, except the `SEA`, must have a way for
    ## water to leave it.
    hydro_nodes = _hydro_nodes(data)
    has_outflow =
        union(data.hydro_stations.HEAD_WATER_FROM, data.hydro_arcs.ORIG)
    for node in hydro_nodes
        @assert node == "SEA" || node in has_outflow "No outflow from $node"
    end
    @assert data.reservoirs.RESERVOIR ⊆ hydro_nodes
    @assert data.inflows.CATCHMENT ⊆ hydro_nodes
    ## Every catchment must have an inflow in every week of every historical
    ## year.
    @assert size(data.inflows, 1) ==
            length(unique(data.inflows.CATCHMENT)) *
            length(unique(data.inflows.YEAR)) *
            length(weeks)
    ## The historical years must be consecutive, so that the historical
    ## sampling schemes can step from one year to the next.
    years = sort(unique(data.inflows.YEAR))
    @assert years == years[1]:years[end] "Historical years are not consecutive"
    ## Water must be able to reach the `SEA` from every reservoir and every
    ## catchment with an inflow using only hydro arcs. Because the bounds on the
    ## hydro arcs are soft constraints (see `_add_hydro_network`), this ensures
    ## that every subproblem is feasible: water that cannot be stored or sent
    ## through a turbine can always be spilled.
    drains_to_sea = Set(["SEA"])
    while true
        new_nodes = setdiff(
            [
                r.ORIG for
                r in eachrow(data.hydro_arcs) if r.DEST in drains_to_sea
            ],
            drains_to_sea,
        )
        if isempty(new_nodes)
            break
        end
        union!(drains_to_sea, new_nodes)
    end
    for node in union(data.reservoirs.RESERVOIR, data.inflows.CATCHMENT)
        @assert node in drains_to_sea "$node cannot spill to the SEA"
    end
    return
end

# ### Helper functions
#
# The following helper functions are used by both the validation and the model.
#
# The nodes of the hydro network are every node that is the end of a hydro arc
# or a hydro station:

function _hydro_nodes(data)
    return unique(
        vcat(
            data.hydro_stations.HEAD_WATER_FROM,
            data.hydro_stations.TAIL_WATER_TO,
            data.hydro_arcs.ORIG,
            data.hydro_arcs.DEST,
        ),
    )
end

# A transmission line is a directed `(from, to)` pair of nodes. Each line has
# one row in the `transmission` table for each of its loss tranches:

function _transmission_lines(data)
    return unique((r.FROM_NODE, r.TO_NODE) for r in eachrow(data.transmission))
end

# `_transmission_losses` returns a vector with a `(capacity, loss_fraction)`
# tuple for each tranche of `line`, in order of the tranche number:
function _transmission_losses(data, line::Tuple{String,String})
    rows = filter(r -> (r.FROM_NODE, r.TO_NODE) == line, data.transmission)
    @assert sort(rows.TRANCHE) == 1:size(rows, 1) "Tranches of $line"
    sort!(rows, :TRANCHE)
    return collect(zip(rows.CAPACITY, rows.LOSS_FRACTION))
end

# ## Configuration
#
# The options that a user can change are collected in `Config`. Everything else
# is a fixed part of the model: to change it, edit this file.
#
# Three of the options are penalties and scalings that are not part of the
# data. Some of the constraints in the hydro network, such as minimum
# environmental flows, can be impossible to satisfy in a drought. To ensure
# that every subproblem is feasible, the bounds in the `hydro_arcs` table are
# soft constraints. A violation is penalized as if it were electricity at
# `penalty_below_min_flow` or `penalty_above_max_flow` (in \$/MWh), where the
# conversion from flow to electricity uses the largest specific power of any
# water in the system (see `_max_specific_power`). The defaults are the
# penalties used by JADE.
#
# The objective is measured in millions of dollars, `objective_scale = 1e-6`,
# to improve the numerical stability of the subproblems.
#
# The uncertainty in each week is the inflows from that week of the historical
# years. By default, every year in the `inflows` table is a realization. Pass
# `inflow_years` to train the policy on a subset of the years, for example, to
# hold some years out for an out-of-sample simulation.

"""
    Config(;
        optimizer,
        current_week::Int = 1,
        discount_factor::Float64 = 0.95,
        penalty_below_min_flow::Float64 = 500.0,
        penalty_above_max_flow::Float64 = 50.0,
        objective_scale::Float64 = 1e-6,
        inflow_years::Union{Nothing,Vector{Int}} = nothing,
    )

The configuration options for `SDDP.PolicyGraph(data, config)`. The keyword
arguments are:

 * `optimizer`: the JuMP optimizer used to solve the subproblems, for example,
   `HiGHS.Optimizer`.
 * `current_week`: the week of the year, from `1` to `52`, in which the model
   starts. The reservoirs start at `INITIAL_LEVEL` in this week.
 * `discount_factor`: the annual discount factor, applied each time the policy
   completes a year.
 * `penalty_below_min_flow`: the cost (\$/MWh) of flow below the minimum flow of
   a hydro arc.
 * `penalty_above_max_flow`: the cost (\$/MWh) of flow above the maximum flow of
   a hydro arc.
 * `objective_scale`: the stage objective, in dollars, is multiplied by this
   before it is passed to the solver.
 * `inflow_years`: the historical years whose inflows are the in-sample
   realizations of the uncertainty in each week. If `nothing`, every year in
   the data is used. Otherwise, a year is used if and only if it is in
   `inflow_years`, and every year in `inflow_years` must be in the data. Any
   `AbstractVector{Int}`, such as `2000:2008`, is converted to a `Vector{Int}`.
"""
Base.@kwdef struct Config
    optimizer::Any
    current_week::Int = 1
    discount_factor::Float64 = 0.95
    penalty_below_min_flow::Float64 = 500.0
    penalty_above_max_flow::Float64 = 50.0
    objective_scale::Float64 = 1e-6
    inflow_years::Union{Nothing,Vector{Int}} = nothing
end

# `inflow_years` can be checked only once we have the data:

function _validate_inflow_years(data::Data, config::Config)
    if config.inflow_years === nothing
        return
    end
    @assert !isempty(config.inflow_years) "`inflow_years` is empty"
    missing_years = setdiff(config.inflow_years, data.inflows.YEAR)
    @assert isempty(missing_years) "No inflows for years: $missing_years"
    return
end

# ## The policy graph
#
# The nodes in the policy graph _always_ correspond to the week of the year:
# node `t` is week `t`, and node `1` is the first week of January. We support
# only cyclic graphs, so the model is infinite horizon. We don't support
# multi-year graphs.
#
# The root node `0` connects to `start_week`. Each week connects to the next,
# and week 52 connects to week 1. The arc that completes a year (the arc into
# `start_week`) has probability `discount_factor`. This means that the policy
# has a `1 - discount_factor` chance of terminating each year, which is
# equivalent to discounting future costs by `discount_factor` each year.

function _create_graph(; start_week::Int, discount_factor::Float64)
    T::Int = 52
    @assert 1 <= start_week <= T
    @assert 0.0 < discount_factor < 1.0
    graph = SDDP.Graph(0)
    for t in 1:T
        SDDP.add_node(graph, t)
    end
    SDDP.add_edge(graph, 0 => start_week, 1.0)
    for t in 1:T
        p = t == start_week ? discount_factor : 1.0
        SDDP.add_edge(graph, mod1(t - 1, T) => t, p)
    end
    return graph
end

# ## Unit conversions
#
# The water balance is in million cubic metres, but flows are in cumecs. A flow
# of one cumec for one hour is 3600 cubic metres.

const MM3_PER_CUMEC_HOUR = 3600 / 1e6

# The specific power of a node is the power (MW) produced by one cumec of water
# as it flows from that node to the `SEA`, taking the path that produces the
# most power. We use it in two places: to convert the flow penalties in
# `Config` from \$/MWh to \$/cumec, using the largest specific power of any
# water in the system, so that the penalty for violating a bound by one cumec
# is the same at every arc in the network; and to convert the marginal value
# of water in a reservoir from \$/Mm³ to \$/MWh.

function _specific_power(data::Data)
    specific_power = Dict{String,Float64}("SEA" => 0.0)
    function _specific_power(node)
        if !haskey(specific_power, node)
            value = 0.0
            for r in eachrow(data.hydro_stations)
                if r.HEAD_WATER_FROM == node
                    value = max(
                        value,
                        r.SPECIFIC_POWER + _specific_power(r.TAIL_WATER_TO),
                    )
                end
            end
            for r in eachrow(data.hydro_arcs)
                if r.ORIG == node
                    value = max(value, _specific_power(r.DEST))
                end
            end
            specific_power[node] = value
        end
        return specific_power[node]
    end
    foreach(_specific_power, _hydro_nodes(data))
    return specific_power
end

_max_specific_power(data::Data) = maximum(values(_specific_power(data)))

# ## Creating the model
#
# The main entry point of JADE is a new method of `SDDP.PolicyGraph`
# that takes our `Data` and `Config`. This is not type piracy, because we own
# the types of both arguments.
#
# It builds the policy graph, and then builds the subproblem in each node in
# three steps: the hydro network, the power network, and the objective. It then
# stores the factors that `simulate` needs to convert the duals of the
# subproblem into prices (see "Recording prices and water values" below). The
# networks are described in the next section. Each adds its own cost
# expression to the subproblem, so the objective is just their sum.

"""
    SDDP.PolicyGraph(data::Data, config::Config)

Create an `SDDP.PolicyGraph` of the JADE model from `data` and the options in
`config`.
"""
function SDDP.PolicyGraph(data::Data, config::Config)
    _validate_inflow_years(data, config)
    graph = _create_graph(;
        start_week = config.current_week,
        discount_factor = config.discount_factor,
    )
    return SDDP.PolicyGraph(
        graph;
        sense = :Min,
        lower_bound = 0.0,
        optimizer = config.optimizer,
    ) do sp, week
        _add_hydro_network(sp, data, config, week)
        _add_power_network(sp, data, week)
        SDDP.@stageobjective(
            sp,
            config.objective_scale *
            (sp[:e_penalty_cost] + sp[:e_thermal_cost] + sp[:e_lost_load_cost]),
        )
        _add_dual_conversions(sp, data, config, week)
        return
    end
end

# ## The subproblem
#
# The subproblem in each week is built from two networks:
#
#  1. A network for the flows of water. Each node is a reservoir, which can
#     store water between weeks, or a junction, which cannot store water. Water
#     flows between nodes through the turbines of hydro stations (which produce
#     power) and through hydro arcs (which do not).
#  2. A network for the flows of power. Each node has a demand, which must be
#     met by generation from hydro, thermal, and fixed stations, by transmission
#     from other nodes, or by shedding load.
#
# Each week is divided into load blocks. Decisions are made for each load
# block, with the exception of the reservoir levels, which are balanced only
# over the whole week.
#
# Each network is a function that adds variables, constraints, and a cost
# expression to the JuMP model `sp`. The power network needs the hydro
# generation from the hydro network; it looks it up with `sp[:name]`.
#
# The naming convention is that state variables start with `x_`, control
# variables start with `u_`, random variables start with `w_`, and expressions
# start with `e_`.
#
# ### Helpers
#
# Many of the tables have one row for each week. `_this_week` picks out the
# rows for the current week. Both networks need the set of load blocks `B`,
# the hours in each load block, and the outages in this week. The hydro
# network needs the in-sample inflows of the week: one `Dict` of catchment to
# inflow for each year in `config.inflow_years`.

function _this_week(df::DataFrames.DataFrame, week::Int)
    return filter(r -> r.WEEK == week, df)
end

_load_blocks(data::Data) = unique(data.hours_per_block.LOAD_BLOCK)

function _hours(data::Data, week::Int)
    rows = _this_week(data.hours_per_block, week)
    return Dict(r.LOAD_BLOCK => r.HOURS for r in eachrow(rows))
end

function _outages(data::Data, week::Int)
    rows = _this_week(data.station_outages, week)
    return Dict(r.GENERATOR => r.OUTAGE for r in eachrow(rows))
end

function _inflow_realizations(data::Data, config::Config, week::Int)
    df = _this_week(data.inflows, week)
    if config.inflow_years !== nothing
        filter!(r -> r.YEAR in config.inflow_years, df)
    end
    return [
        Dict(g.CATCHMENT .=> g.INFLOW) for g in DataFrames.groupby(df, :YEAR)
    ]
end

# ### The hydro network
#
# The sets are `N`, the nodes; `R`, the reservoirs; `J`, the junctions; `A`,
# the hydro arcs; and `H`, the hydro stations. The hydro network adds the
# penalty for violating the soft bounds on the hydro arcs, `e_penalty_cost`,
# and the power produced by each hydro station, `e_hydro_generation`.

function _add_hydro_network(
    sp::JuMP.Model,
    data::Data,
    config::Config,
    week::Int,
)
    B, hours, outage =
        _load_blocks(data), _hours(data, week), _outages(data, week)
    N = _hydro_nodes(data)
    R = data.reservoirs.RESERVOIR
    J = setdiff(N, R, ["SEA"])
    A = [(r.ORIG, r.DEST) for r in eachrow(data.hydro_arcs)]
    H = data.hydro_stations.GENERATOR
    reservoir = Dict(r.RESERVOIR => r for r in eachrow(data.reservoirs))
    arc = Dict((r.ORIG, r.DEST) => r for r in eachrow(data.hydro_arcs))
    hydro = Dict(r.GENERATOR => r for r in eachrow(data.hydro_stations))
    ## The storage in reservoir `r` at the start and end of the week.
    ## Unit: million cubic metres [Mm³]
    JuMP.@variable(
        sp,
        0 <= x_storage[r in R] <= reservoir[r].MAX_LEVEL,
        SDDP.State,
        initial_value = reservoir[r].INITIAL_LEVEL,
    )
    ## The inflow into node `n`. Nodes without an inflow in the data have zero
    ## inflow. Each realization of the uncertainty is the inflows from this
    ## week of one historical year in `config.inflow_years`.
    ## Unit: cumecs [m³/s]
    JuMP.@variable(sp, w_inflow[n in N] == 0)
    Ω = _inflow_realizations(data, config, week)
    SDDP.parameterize(sp, Ω) do ω
        for (n, inflow) in ω
            JuMP.fix(w_inflow[n], inflow)
        end
        return
    end
    ## The flow along hydro arc `a`, and the flow through the turbine of hydro
    ## station `h`, in load block `b`.
    ## Unit: cumecs [m³/s]
    JuMP.@variable(sp, u_flow[a in A, b in B] >= 0)
    JuMP.@variable(sp, u_turbine[h in H, b in B] >= 0)
    ## The bounds on the hydro arcs are soft constraints. `u_flow_below` is the
    ## violation of the minimum flow, and `u_flow_above` is the violation of the
    ## maximum flow.
    ## Unit: cumecs [m³/s]
    JuMP.@variable(sp, u_flow_below[a in A, b in B] >= 0)
    JuMP.@variable(sp, u_flow_above[a in A, b in B] >= 0)
    JuMP.@constraint(
        sp,
        [a in A, b in B],
        u_flow[a, b] + u_flow_below[a, b] >= arc[a].MIN_FLOW,
    )
    JuMP.@constraint(
        sp,
        [a in A, b in B; !ismissing(arc[a].MAX_FLOW)],
        u_flow[a, b] - u_flow_above[a, b] <= arc[a].MAX_FLOW,
    )
    ## The net flow into node `n` in load block `b`: the inflow, plus the flow
    ## in from arcs and turbines, minus the flow out to arcs and turbines.
    ## Unit: cumecs [m³/s]
    JuMP.@expression(
        sp,
        e_net_flow[n in N, b in B],
        w_inflow[n] + sum(u_flow[(i, j), b] for (i, j) in A if j == n) -
        sum(u_flow[(i, j), b] for (i, j) in A if i == n) +
        sum(u_turbine[h, b] for h in H if hydro[h].TAIL_WATER_TO == n) -
        sum(u_turbine[h, b] for h in H if hydro[h].HEAD_WATER_FROM == n)
    )
    ## Reservoirs can store water. We don't try to balance the reservoirs
    ## within a week; we ensure only that they are within bounds at the start
    ## and end of the week.
    JuMP.@constraint(
        sp,
        c_reservoir_balance[r in R],
        x_storage[r].out - x_storage[r].in ==
        sum(MM3_PER_CUMEC_HOUR * hours[b] * e_net_flow[r, b] for b in B),
    )
    ## Junctions cannot store water, so flow in equals flow out in every load
    ## block. The `SEA` is the sink of the network, so it has no balance.
    JuMP.@constraint(
        sp,
        c_junction_balance[j in J, b in B],
        e_net_flow[j, b] == 0
    )
    ## The power produced by a hydro station is linear in the flow through its
    ## turbine, and it is limited by the capacity of the station, less any
    ## outage.
    ## Unit: MW
    JuMP.@expression(
        sp,
        e_hydro_generation[h in H, b in B],
        hydro[h].SPECIFIC_POWER * u_turbine[h, b],
    )
    JuMP.@constraint(
        sp,
        [h in H, b in B],
        e_hydro_generation[h, b] <= hydro[h].CAPACITY - get(outage, h, 0.0),
    )
    ## The cost of violating the soft bounds, converted from $/MWh to $/cumec
    ## with the largest specific power in the system.
    ## Unit: $
    max_specific_power = _max_specific_power(data)
    JuMP.@expression(
        sp,
        e_penalty_cost,
        sum(
            hours[b] *
            max_specific_power *
            (
                config.penalty_below_min_flow * u_flow_below[a, b] +
                config.penalty_above_max_flow * u_flow_above[a, b]
            ) for a in A, b in B
        ),
    )
    return
end

# ### The power network
#
# The sets are `P`, the nodes; `G`, the thermal stations; `T`, the tranches of
# the transmission lines; and `K`, the tranches of lost load. The power network
# adds the cost of thermal generation, `e_thermal_cost`, and of shedding load,
# `e_lost_load_cost`.

function _add_power_network(sp::JuMP.Model, data::Data, week::Int)
    B, hours, outage =
        _load_blocks(data), _hours(data, week), _outages(data, week)
    P = unique(data.demand.NODE)
    G = data.thermal_stations.GENERATOR
    T = [
        (r.FROM_NODE, r.TO_NODE, r.TRANCHE) for r in eachrow(data.transmission)
    ]
    K = [(r.NODE, r.COST) for r in eachrow(data.lost_load)]
    thermal = Dict(r.GENERATOR => r for r in eachrow(data.thermal_stations))
    tranche = Dict(
        (r.FROM_NODE, r.TO_NODE, r.TRANCHE) => r for
        r in eachrow(data.transmission)
    )
    lost_load = Dict((r.NODE, r.COST) => r for r in eachrow(data.lost_load))
    demand = Dict(
        (r.NODE, r.LOAD_BLOCK) => r.DEMAND for
        r in eachrow(_this_week(data.demand, week))
    )
    ## The output of the fixed stations at each node, which is the same in
    ## every week.
    ## Unit: MW
    fixed = DataFrames.innerjoin(
        data.fixed_station_output,
        data.fixed_stations;
        on = :STATION,
    )
    fixed_generation = Dict(
        (p, b) => sum(
            r.OUTPUT for
            r in eachrow(fixed) if r.NODE == p && r.LOAD_BLOCK == b;
            init = 0.0,
        ) for p in P, b in B
    )
    ## The hydro generation, from the hydro network.
    hydro_node =
        Dict(r.GENERATOR => r.NODE for r in eachrow(data.hydro_stations))
    e_hydro_generation = sp[:e_hydro_generation]
    ## The generation of thermal station `g` in load block `b`, limited by the
    ## capacity of the station, less any outage.
    ## Unit: MW
    JuMP.@variable(
        sp,
        0 <=
        u_thermal[g in G, b in B] <=
        thermal[g].CAPACITY - get(outage, g, 0.0),
    )
    ## The flow in tranche `t` of a transmission line in load block `b`. Because
    ## the loss fractions of the tranches of a line are increasing, the cheaper
    ## tranches fill first and the losses are piecewise linear.
    ## Unit: MW
    JuMP.@variable(
        sp,
        0 <= u_transmission[t in T, b in B] <= tranche[t].CAPACITY
    )
    ## The demand shed in lost load tranche `k` in load block `b`.
    ## Unit: MW
    JuMP.@variable(
        sp,
        0 <=
        u_lost_load[k in K, b in B] <=
        lost_load[k].PROPORTION * demand[(k[1], b)],
    )
    ## Supply must meet demand at each node in each load block. Half of the
    ## losses on a transmission line are attributed to each end of the line.
    ## The constraint is `>=` rather than `==` so that excess generation, for
    ## example, from minimum flows through a turbine, can be freely disposed of.
    ## Unit: MW
    JuMP.@constraint(
        sp,
        c_demand[p in P, b in B],
        sum(
            e_hydro_generation[h, b] for
            h in keys(hydro_node) if hydro_node[h] == p
        ) +
        sum(u_thermal[g, b] for g in G if thermal[g].NODE == p) +
        fixed_generation[(p, b)] +
        sum(
            (1 - 0.5 * tranche[t].LOSS_FRACTION) * u_transmission[t, b] for
            t in T if t[2] == p
        ) - sum(
            (1 + 0.5 * tranche[t].LOSS_FRACTION) * u_transmission[t, b] for
            t in T if t[1] == p
        ) + sum(u_lost_load[k, b] for k in K if k[1] == p) >= demand[(p, b)]
    )
    ## The cost of thermal generation is the operations and maintenance cost,
    ## plus the cost of the fuel and its carbon emissions.
    ## Unit: $/MWh
    co2_content = Dict(r.FUEL => r.CO2_CONTENT for r in eachrow(data.fuels))
    co2_price = only(_this_week(data.co2_prices, week).CO2_PRICE)
    fuel_cost = Dict(
        r.FUEL => r.PRICE + co2_content[r.FUEL] * co2_price for
        r in eachrow(_this_week(data.fuel_prices, week))
    )
    thermal_cost = Dict(
        g =>
            thermal[g].OMCOST +
            thermal[g].HEAT_RATE * fuel_cost[thermal[g].FUEL] for g in G
    )
    ## Unit: $
    JuMP.@expression(
        sp,
        e_thermal_cost,
        sum(hours[b] * thermal_cost[g] * u_thermal[g, b] for g in G, b in B),
    )
    JuMP.@expression(
        sp,
        e_lost_load_cost,
        sum(
            hours[b] * lost_load[k].COST * u_lost_load[k, b] for k in K, b in B
        ),
    )
    return
end

# ## Training the policy
#
# The policy graph is cyclic, so a forward pass of `SDDP.train` continues until
# the discounting terminates it, which takes 20 years on average. In the early
# iterations, the cost-to-go functions are poor, so a long forward pass adds
# cuts at states that the policy is unlikely to visit again. Instead, the
# default sampling scheme limits the forward pass to one year (52 weeks) for
# the first 104 iterations, then to two years for the next 104 iterations, and
# so on. A forward pass can still end earlier, each time it completes a year,
# because of the discounting.

"""
    default_sampling_scheme()

Return the default sampling scheme of [`train`](@ref): an
`SDDP.InSampleMonteCarlo` in which the forward pass of iteration `i` is limited
to `52 * cld(i, 104)` weeks.

The iteration count starts at `1` each time this function is called, so pass
the same sampling scheme to every call of `train` that continues training the
same model.
"""
function default_sampling_scheme()
    return SDDP.InSampleMonteCarlo(; rollout_limit = i -> 52 * cld(i, 104))
end

"""
    train(
        model::SDDP.PolicyGraph;
        sampling_scheme::SDDP.AbstractSamplingScheme = default_sampling_scheme(),
        kwargs...,
    )

Train the policy in `model`, which must have been created by
`SDDP.PolicyGraph(data, config)`, with `SDDP.train`.

The default `sampling_scheme` is [`default_sampling_scheme`](@ref). The other
keyword arguments are passed to `SDDP.train`, for example,
`iteration_limit = 100`.
"""
function train(
    model::SDDP.PolicyGraph;
    sampling_scheme::SDDP.AbstractSamplingScheme = default_sampling_scheme(),
    kwargs...,
)
    return SDDP.train(model; sampling_scheme, kwargs...)
end

# ## Simulating the policy
#
# Once the policy is trained, `SDDP.simulate` evaluates it. This section adds
# two utilities: sampling schemes that replay the historical inflows, and
# `simulate`, a wrapper around `SDDP.simulate` that records prices and water
# values.
#
# ### The historical inflows
#
# During training, the inflows in each week are sampled independently from the
# historical years. To evaluate a policy, we often want to simulate it on the
# historical sequence of inflows instead. There are two ways to do so, and both
# return an `SDDP.Historical` sampling scheme, which replays the scenarios it is
# given in order.
#
# Every scenario starts in week `config.current_week`, because that is the first
# week of the model, and runs for a whole number of years. A scenario that
# starts after week 1 ends in week `config.current_week - 1` of a later year.
# The scenarios never wrap around from the last historical year to the first, so
# if the model starts in week 39, weeks 1 to 38 of the first historical year and
# weeks 39 to 52 of the last are never used.
#
# Both samplers look up the inflows of a week by `(year, week)`:

function _inflows_by_week(data::Data)
    groups = DataFrames.groupby(data.inflows, [:YEAR, :WEEK])
    return Dict(
        (key.YEAR, key.WEEK) => Dict(df.CATCHMENT .=> df.INFLOW) for
        (key, df) in pairs(groups)
    )
end

# and build a scenario as a vector of `(node, noise)` tuples, where the node of
# the policy graph is the week:

function _historical_scenario(
    inflows::Dict,
    start_year::Int,
    start_week::Int,
    number_weeks::Int,
)
    return map(0:(number_weeks-1)) do i
        year = start_year + fld(start_week - 1 + i, 52)
        week = mod1(start_week + i, 52)
        return (week, inflows[(year, week)])
    end
end

"""
    HistoricalYears(data::Data, config::Config)

Return an `SDDP.Historical` sampling scheme with one scenario of 52 weeks for
each historical year in `data`, for simulating a model created by
`SDDP.PolicyGraph(data, config)`.

Each scenario starts in week `config.current_week` of a historical year. If the
model starts after week 1, the scenario ends in the following year, so there is
one fewer scenario than there are historical years.

The scenarios are sampled in order, so simulate one replication for each:

```julia
sampling_scheme = JADE.HistoricalYears(data, config)
SDDP.simulate(model, length(sampling_scheme.scenarios); sampling_scheme)
```
"""
function HistoricalYears(data::Data, config::Config)
    inflows = _inflows_by_week(data)
    years = sort(unique(data.inflows.YEAR))
    ## A scenario that starts after week 1 needs the following year, so the
    ## last historical year cannot start a scenario.
    last_year = config.current_week == 1 ? years[end] : years[end] - 1
    scenarios = map(years[1]:last_year) do year
        return _historical_scenario(inflows, year, config.current_week, 52)
    end
    return SDDP.Historical(scenarios)
end

"""
    HistoricalSequence(
        data::Data,
        config::Config;
        start_year::Int = minimum(data.inflows.YEAR),
        end_year::Int = maximum(data.inflows.YEAR),
    )

Return an `SDDP.Historical` sampling scheme with a single scenario that starts
in week `config.current_week` of `start_year` and runs for a whole number of
years, using the historical inflows in `data` up to and including `end_year`,
for simulating a model created by `SDDP.PolicyGraph(data, config)`.

If the model starts in week 1, the scenario ends in week 52 of `end_year`.
Otherwise, it ends in week `config.current_week - 1` of `end_year`. The
defaults use the whole historical record.

```julia
config = JADE.Config(; optimizer = HiGHS.Optimizer, current_week = 39)
sampling_scheme = JADE.HistoricalSequence(
    data,
    config;
    start_year = 1980,
    end_year = 1990,
)
SDDP.simulate(model, 1; sampling_scheme)
```
"""
function HistoricalSequence(
    data::Data,
    config::Config;
    start_year::Int = minimum(data.inflows.YEAR),
    end_year::Int = maximum(data.inflows.YEAR),
)
    inflows = _inflows_by_week(data)
    years = sort(unique(data.inflows.YEAR))
    @assert start_year in years "No historical inflows in $start_year"
    @assert end_year in years "No historical inflows in $end_year"
    ## The number of whole years from `start_year` to the end of `end_year`.
    number_years = end_year - start_year + (config.current_week == 1)
    @assert(
        number_years >= 1,
        "No whole years between $start_year and $end_year",
    )
    scenario = _historical_scenario(
        inflows,
        start_year,
        config.current_week,
        52 * number_years,
    )
    return SDDP.Historical([scenario])
end

# ### Recording prices and water values
#
# `SDDP.simulate` records the value of any variable. In a hydro-thermal model,
# two duals are just as interesting: the price of energy at each node, which is
# the dual of the demand constraint, and the marginal value of water in each
# reservoir, which is the dual of the reservoir balance. Both need their units
# converted, because the objective is scaled by `objective_scale`, the demand
# constraint is in MW for a whole load block, and the reservoir balance is in
# Mm³.
#
# The conversion factors depend on the data and the configuration, but we do
# not want `simulate` to take either as an argument, because then a model could
# be simulated with data it was not built from. Instead, `SDDP.PolicyGraph`
# stores the factors in the `ext` dictionary of each subproblem when it builds
# the model:

function _add_dual_conversions(
    sp::JuMP.Model,
    data::Data,
    config::Config,
    week::Int,
)
    ## The dual of `c_demand[p, b]` is the scaled cost of one more MW of demand
    ## for the whole of load block `b`. Dividing by the hours in the block gives
    ## $/MWh.
    hours = _hours(data, week)
    sp.ext[:price_scale] =
        Dict(b => 1 / (config.objective_scale * h) for (b, h) in hours)
    ## The dual of `c_reservoir_balance[r]` is the scaled change in cost from
    ## one more Mm³ of water in reservoir `r`. It is negative, because more
    ## water reduces the cost. One Mm³ of water in `r` can produce
    ## `specific_power[r] / MM3_PER_CUMEC_HOUR` MWh.
    specific_power = _specific_power(data)
    scale = -MM3_PER_CUMEC_HOUR / config.objective_scale
    sp.ext[:water_value_scale] =
        Dict(r => scale / specific_power[r] for r in data.reservoirs.RESERVOIR)
    return
end

# The two recorders read the duals from the subproblem and apply the factors:

function _price(sp::JuMP.Model)
    price = JuMP.dual.(sp[:c_demand])
    for p in axes(price, 1), b in axes(price, 2)
        price[p, b] *= sp.ext[:price_scale][b]
    end
    return price
end

function _marginal_water_value(sp::JuMP.Model)
    water_value = JuMP.dual.(sp[:c_reservoir_balance])
    for r in axes(water_value, 1)
        water_value[r] *= sp.ext[:water_value_scale][r]
    end
    return water_value
end

"""
    simulate(model::SDDP.PolicyGraph, number_replications::Int; kwargs...)

Simulate the policy in `model`, which must have been created by
`SDDP.PolicyGraph(data, config)`, for `number_replications` replications, and
return the vector of simulations from `SDDP.simulate`.

In addition to what `SDDP.simulate` always records, each stage of each
simulation contains:

 * `:x_storage`: the storage in each reservoir at the start (`.in`) and the end
   (`.out`) of the week (Mm³)
 * `:u_thermal`: the generation of each thermal station in each load block (MW)
 * `:price`: the price of energy at each node in each load block (\$/MWh), which
   is the dual of the demand constraint
 * `:marginal_water_value`: the value of one more MWh of water in each
   reservoir (\$/MWh), which is the dual of the reservoir balance constraint,
   converted from Mm³ to MWh with the specific power of the reservoir

The keyword arguments are passed to `SDDP.simulate`, for example,
`sampling_scheme = HistoricalYears(data, config)`.
"""
function simulate(model::SDDP.PolicyGraph, number_replications::Int; kwargs...)
    return SDDP.simulate(
        model,
        number_replications,
        [:x_storage, :u_thermal];
        custom_recorders = Dict{Symbol,Function}(
            :price => _price,
            :marginal_water_value => _marginal_water_value,
        ),
        kwargs...,
    )
end

# ## Converting CSV files to a database
#
# CSV files are a convenient format for editing the input data by hand. This
# section converts a directory of CSV files into a database with the schema
# above. Each CSV file becomes one table. Most files are a direct copy
# of the table, but some are in a "wide" format that is easier to edit, for
# example, with one column per load block.
#
# Missing values in the CSV files are written as `NA`.
#
# `data` is the CSV source of `data.db`. If you change the CSV files, regenerate
# the database: the tests check that the two are the same.
#
# ### Converting a directory
#
# The keyword constructor `Data(; filename, directory, force)` converts the CSV
# files in `directory` into the database `filename`, and then reads it back
# with `Data(filename)`. Reading it back is what validates the data, so a
# database that fails validation is deleted rather than left on disk.

function Data(;
    filename::String,
    directory::Union{Nothing,String} = nothing,
    force::Bool = false,
)
    if directory === nothing
        return Data(filename)
    end
    if isfile(filename)
        @assert force "$filename exists. Use `force = true` to overwrite it"
        rm(filename)
    end
    try
        _csv_to_db(directory, filename)
        return Data(filename)
    catch
        rm(filename; force = true)  # `force` here means "if it exists"
        rethrow()
    end
end

# `_csv_to_db` reads each CSV file with its own reader (see below) and writes
# it to the table of the same name. The readers are in the order of `TABLES`.

function _csv_to_db(directory::String, filename::String)
    readers = [
        _csv_hours_per_block,
        _csv_demand,
        _csv_lost_load,
        _csv_fixed_stations,
        _csv_fixed_station_output,
        _csv_thermal_stations,
        _csv_fuels,
        _csv_fuel_prices,
        _csv_co2_prices,
        _csv_hydro_stations,
        _csv_hydro_arcs,
        _csv_reservoirs,
        _csv_inflows,
        _csv_station_outages,
        _csv_transmission,
    ]
    db = SQLite.DB(filename)
    try
        _create_tables(db)
        for (table, reader) in zip(TABLES, readers)
            SQLite.load!(reader(directory), db, table)
        end
        ## Rebuild the file to remove the free space left by inserting rows that
        ## are not in the order of the primary key.
        DBInterface.execute(db, "VACUUM")
    finally
        DBInterface.close!(db)
    end
    return
end

# ### Reading a CSV file
#
# `_read_csv` reads `directory/filename`. If `columns` is provided, the file
# must have exactly these columns.
function _read_csv(directory::String, filename::String, columns = nothing)
    path = joinpath(directory, filename)
    @assert isfile(path) "Missing input file: $path"
    df = CSV.read(
        path,
        DataFrames.DataFrame;
        missingstring = "NA",
        stringtype = String,
    )
    if columns !== nothing
        @assert names(df) == columns "The columns of $filename must be $columns"
    end
    return df
end

# `_stack` converts a wide `df`, which has the `id_columns` followed by one
# column for each value of `variable_name`, into the long format of the
# database.
function _stack(df, id_columns; variable_name, value_name)
    @assert(
        names(df)[1:length(id_columns)] == id_columns,
        "The first columns must be $id_columns",
    )
    @assert size(df, 2) > length(id_columns)
    return DataFrames.stack(
        df,
        DataFrames.Not(id_columns);
        variable_name,
        value_name,
    )
end

# ### The CSV files
#
# There is one reader for each table in the schema, in the same order as the
# schema. Each reader returns a `DataFrames.DataFrame` with the same columns as
# its table, so `_csv_to_db` can write it directly.
#
# Most of the files are a copy of their table: they must have exactly the
# columns of the table, in the same order, and the units and the allowed values
# of each column are those given in the schema. Names, such as nodes, stations,
# and fuels, are case-sensitive and must match across files.
#
# The other files are in a wide format, with one column for each load block or
# each catchment. The names of those columns are not fixed: they become the
# values of the `LOAD_BLOCK` or `CATCHMENT` column of the table. Every wide
# file with load blocks must use the same set of load blocks.
#
# The code block under each file is an example of the structure of the file:
# the header row, followed by the first few rows of the default data. The real
# files have many more rows, and the examples are not valid on their own.

# #### `hours_per_block.csv`
#
# Wide. One row for each `WEEK` from 1 to 52. The first column is `WEEK`, and
# every other column is a load block. The values are the hours in that load
# block of that week; they should sum to the hours in the week.
#
# ```
# WEEK,PEAK,SHOULDER,OFFPEAK
# 1,35,76,57
# 2,59.5,56,52.5
# ```

function _csv_hours_per_block(directory::String)
    df = _read_csv(directory, "hours_per_block.csv")
    return _stack(
        df,
        ["WEEK"];
        variable_name = :LOAD_BLOCK,
        value_name = :HOURS,
    )
end

# #### `demand.csv`
#
# Wide. One row for each `NODE` of the power network and each `WEEK`. The first
# two columns are `NODE` and `WEEK`, and every other column is a load block.
# The values are the average demand (MW) in that load block. Every node of the
# power network must appear here, even if its demand is zero: this file defines
# the set of nodes.
#
# ```
# NODE,WEEK,PEAK,SHOULDER,OFFPEAK
# HAY,1,293.8821429,161.5040789,114.590614
# HAY,2,276.5505882,201.8775,131.1262857
# ```

function _csv_demand(directory::String)
    df = _read_csv(directory, "demand.csv")
    return _stack(
        df,
        ["NODE", "WEEK"];
        variable_name = :LOAD_BLOCK,
        value_name = :DEMAND,
    )
end

# #### `lost_load.csv`
#
# A copy of the `lost_load` table: `NODE`, `COST` (\$/MWh), and `PROPORTION`
# of the demand that can be shed at that cost. One row for each tranche. A node
# may have several tranches at different costs, and their proportions must sum
# to at least 1, so that all demand can be shed.
#
# ```
# NODE,COST,PROPORTION
# NI,1000,1.0
# HAY,500,0.1
# HAY,1000,0.9
# ```

function _csv_lost_load(directory::String)
    columns = ["NODE", "COST", "PROPORTION"]
    return _read_csv(directory, "lost_load.csv", columns)
end

# #### `fixed_stations.csv`
#
# A copy of the `fixed_stations` table: `STATION` and the `NODE` of the power
# network that it generates at.
#
# ```
# STATION,NODE
# WestWind,HAY
# Tararua,NI
# ```

function _csv_fixed_stations(directory::String)
    return _read_csv(directory, "fixed_stations.csv", ["STATION", "NODE"])
end

# #### `fixed_station_output.csv`
#
# Wide. One row for each `STATION` in `fixed_stations.csv`. The first column is
# `STATION`, and every other column is a load block. The values are the output
# (MW) in that load block, which is the same in every week.
#
# ```
# STATION,PEAK,SHOULDER,OFFPEAK
# WestWind,60,60,60
# Tararua,37,37,37
# ```

function _csv_fixed_station_output(directory::String)
    df = _read_csv(directory, "fixed_station_output.csv")
    return _stack(
        df,
        ["STATION"];
        variable_name = :LOAD_BLOCK,
        value_name = :OUTPUT,
    )
end

# #### `thermal_stations.csv`
#
# A copy of the `thermal_stations` table: `GENERATOR`, `NODE`, `FUEL` (a name
# in `fuels.csv`), `HEAT_RATE` (GJ/MWh), `CAPACITY` (MW), and `OMCOST`, the
# operations and maintenance cost (\$/MWh).
#
# ```
# GENERATOR,NODE,FUEL,HEAT_RATE,CAPACITY,OMCOST
# Stratford_220KV,NI,gas,7.6,377,0
# Huntly_e3p,NI,gas,7.2,403,0
# ```

function _csv_thermal_stations(directory::String)
    columns = ["GENERATOR", "NODE", "FUEL", "HEAT_RATE", "CAPACITY", "OMCOST"]
    return _read_csv(directory, "thermal_stations.csv", columns)
end

# #### `fuels.csv`
#
# A copy of the `fuels` table: `FUEL` and its `CO2_CONTENT` (tonnes of CO₂ per
# GJ). One row for each fuel used by a thermal station.
#
# ```
# FUEL,CO2_CONTENT
# coal,0.0912
# gas,0.0528
# ```

function _csv_fuels(directory::String)
    return _read_csv(directory, "fuels.csv", ["FUEL", "CO2_CONTENT"])
end

# #### `fuel_prices.csv`
#
# A copy of the `fuel_prices` table: `WEEK`, `FUEL`, and `PRICE` (\$/GJ). One
# row for each fuel in `fuels.csv` and each week from 1 to 52.
#
# ```
# WEEK,FUEL,PRICE
# 1,coal,5.6
# 1,gas,6.68
# ```

function _csv_fuel_prices(directory::String)
    return _read_csv(directory, "fuel_prices.csv", ["WEEK", "FUEL", "PRICE"])
end

# #### `co2_prices.csv`
#
# A copy of the `co2_prices` table: `WEEK` and `CO2_PRICE` (\$/tonne of CO₂).
# One row for each week from 1 to 52.
#
# ```
# WEEK,CO2_PRICE
# 1,25.0
# 2,25.0
# ```

function _csv_co2_prices(directory::String)
    return _read_csv(directory, "co2_prices.csv", ["WEEK", "CO2_PRICE"])
end

# #### `hydro_stations.csv`
#
# A copy of the `hydro_stations` table: `GENERATOR`; the hydro nodes that the
# turbine takes water from, `HEAD_WATER_FROM`, and returns it to,
# `TAIL_WATER_TO`; the `NODE` of the power network it generates at; its
# `CAPACITY` (MW); and its `SPECIFIC_POWER` (MW per cumec). Spillways are not
# listed here: they are rows in `hydro_arcs.csv`.
#
# ```
# GENERATOR,HEAD_WATER_FROM,TAIL_WATER_TO,NODE,CAPACITY,SPECIFIC_POWER
# Arapuni,Lake_Arapuni,Lake_Karapiro,NI,196.7,0.439847649
# ```

function _csv_hydro_stations(directory::String)
    columns = [
        "GENERATOR",
        "HEAD_WATER_FROM",
        "TAIL_WATER_TO",
        "NODE",
        "CAPACITY",
        "SPECIFIC_POWER",
    ]
    return _read_csv(directory, "hydro_stations.csv", columns)
end

# #### `hydro_arcs.csv`
#
# A copy of the `hydro_arcs` table: the hydro nodes `ORIG` and `DEST`, and the
# `MIN_FLOW` and `MAX_FLOW` (cumec) of the arc. `MAX_FLOW` may be `NA`, which
# means there is no limit. The arcs include rivers, canals, and spillways, and
# together they must let water from every reservoir and every catchment reach
# the `SEA`.
#
# ```
# ORIG,DEST,MIN_FLOW,MAX_FLOW
# Lake_Wanaka,Lake_Dunstan,0.0,NA
# Clyde_tail,Lake_Roxburgh,50,1000
# ```

function _csv_hydro_arcs(directory::String)
    columns = ["ORIG", "DEST", "MIN_FLOW", "MAX_FLOW"]
    return _read_csv(directory, "hydro_arcs.csv", columns)
end

# #### `reservoirs.csv`
#
# A copy of the `reservoirs` table: `RESERVOIR`, a hydro node, and its
# `INITIAL_LEVEL` and `MAX_LEVEL` (Mm³). Hydro nodes that are not listed here
# are junctions, which cannot store water.
#
# ```
# RESERVOIR,INITIAL_LEVEL,MAX_LEVEL
# Lake_Hawea,332.074,1378.764328
# ```

function _csv_reservoirs(directory::String)
    columns = ["RESERVOIR", "INITIAL_LEVEL", "MAX_LEVEL"]
    return _read_csv(directory, "reservoirs.csv", columns)
end

# #### `inflows.csv`
#
# Wide. One row for each historical `YEAR` and `WEEK`, with every year
# complete from week 1 to 52. The first two columns are `YEAR` and `WEEK`, and
# every other column is a catchment, which must be a hydro node. The values are
# the average inflow (cumec) into that catchment in that week.
#
# ```
# YEAR,WEEK,Lake_Aviemore,Lake_Benmore,Lake_Hawea,...
# 1970,1,27,67,71,...
# ```

function _csv_inflows(directory::String)
    df = _read_csv(directory, "inflows.csv")
    df = _stack(
        df,
        ["YEAR", "WEEK"];
        variable_name = :CATCHMENT,
        value_name = :INFLOW,
    )
    return DataFrames.select(df, :CATCHMENT, :YEAR, :WEEK, :INFLOW)
end

# #### `station_outages.csv`
#
# A copy of the `station_outages` table: `GENERATOR`, a hydro or thermal
# station, `WEEK`, and the `OUTAGE` (MW) by which its capacity is reduced.
# List only the weeks with an outage.
#
# ```
# GENERATOR,WEEK,OUTAGE
# Huntly_main_g3,1,250
# Huntly_main_g3,2,250
# ```

function _csv_station_outages(directory::String)
    columns = ["GENERATOR", "WEEK", "OUTAGE"]
    return _read_csv(directory, "station_outages.csv", columns)
end

# #### `transmission.csv`
#
# Wide. One row for each directed line; a line that carries power both ways
# needs two rows. The first three columns are `FROM_NODE`, `TO_NODE`, and
# `CAPACITY` (MW). They are followed by pairs of columns `LOSS_TRANCHE_i` and
# `LOSS_FRACTION_i`, for `i = 1, 2, ...`, which give the piecewise-linear loss
# curve as cumulative breakpoints: the first `LOSS_TRANCHE_1` MW of flow loses
# `LOSS_FRACTION_1` of itself, the flow between `LOSS_TRANCHE_1` and
# `LOSS_TRANCHE_2` loses `LOSS_FRACTION_2`, and so on. The loss fractions must
# increase, and the last `LOSS_TRANCHE_i` must equal the `CAPACITY`. Every row
# has the same number of tranches.
#
# ```
# FROM_NODE,TO_NODE,CAPACITY,LOSS_TRANCHE_1,LOSS_FRACTION_1
# NI,HAY,2343,2343,0.05
# ```
#
# The reader converts the cumulative breakpoints into the width of each tranche,
# which is what the `transmission` table stores.

function _csv_transmission(directory::String)
    df = _read_csv(directory, "transmission.csv")
    @assert names(df)[1:3] == ["FROM_NODE", "TO_NODE", "CAPACITY"]
    n_tranches = div(size(df, 2) - 3, 2)
    @assert n_tranches >= 1 && size(df, 2) == 3 + 2 * n_tranches
    rows = DataFrames.DataFrame(;
        FROM_NODE = String[],
        TO_NODE = String[],
        TRANCHE = Int[],
        CAPACITY = Float64[],
        LOSS_FRACTION = Float64[],
    )
    for r in eachrow(df)
        previous = 0.0
        for i in 1:n_tranches
            limit = r["LOSS_TRANCHE_$i"]
            push!(
                rows,
                (
                    r.FROM_NODE,
                    r.TO_NODE,
                    i,
                    limit - previous,
                    r["LOSS_FRACTION_$i"],
                ),
            )
            previous = limit
        end
        @assert previous == r.CAPACITY "Tranches of $(r.FROM_NODE)-$(r.TO_NODE)"
    end
    return rows
end

end  # module JADE
