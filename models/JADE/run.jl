# Copyright (c) 2017-26: Oscar Dowson and SDDP.jl contributors.          #src
#                                                                        #src
# This Source Code Form is subject to the terms of the Mozilla Public    #src
# License, v2.0. If a copy of the MPL was not distributed with this file #src
# You can obtain one at http://mozilla.org/MPL/2.0/.                     #src

# # Running JADE

# The purpose of this tutorial is to show how to build, train, and simulate
# the [JADE](@ref) model of the New Zealand electricity system.
#
# To keep the tutorial fast enough to run as part of the documentation, we
# train the policy on a small subset of the historical inflows, and we stop
# after a small number of iterations. The settings to change for a real run
# are pointed out as we go.

# ## Packages

# This tutorial requires the following packages:

using SDDP
import HiGHS
import Plots

# ## Loading JADE

# JADE is not part of the `SDDP` module. It is a single file in the SDDP.jl
# repository, which we `include` to define the module `JADE`:

JADE_DIR = joinpath(pkgdir(SDDP), "models", "JADE")
include(joinpath(JADE_DIR, "JADE.jl"));

# ## Data

# The input to JADE is a SQLite database. The repository ships the default
# data as a directory of CSV files, so we convert the directory to a database
# in a temporary location, and then read it back as a `JADE.Data`:

data = JADE.Data(;
    filename = joinpath(mktempdir(), "data.db"),
    directory = joinpath(JADE_DIR, "data"),
)

# Each table in the database is a field of `data`. For example, here are the
# reservoirs, with their initial and maximum levels in million cubic metres:

data.reservoirs

# The historical inflows cover 48 years:

years = sort(unique(data.inflows.YEAR))
(first(years), last(years))

# ## Configuration

# The options of the model are collected in a `JADE.Config`. The only required
# option is the optimizer. We also set `inflow_years` to the last three years of
# the record, so that the uncertainty in each week has three realizations
# instead of 48. This makes the backward pass of SDDP much cheaper. For a real
# run, leave `inflow_years` at its default to use every
# year.

config = JADE.Config(;
    optimizer = HiGHS.Optimizer,
    inflow_years = years[(end-2):end],
)

# ## Training a policy

# `SDDP.PolicyGraph(data, config)` builds the model. The policy graph is cyclic
# with 52 nodes, one for each week of the year:

model = SDDP.PolicyGraph(data, config)

# Here's a visualization of the graph:

## We need `open = false` to build the documentation. Remove if running locally.
SDDP.plot(model, "model_jade.html"; open = false)

# ```@raw html
# <iframe src="../model_jade.html" style="width:100%;height:500px;"></iframe>
# ```

# `JADE.train` wraps `SDDP.train`. Because the graph is cyclic, a forward pass
# continues until the discounting terminates it, which takes 20 years on
# average. To keep the early iterations short, `JADE.train` limits the forward
# pass to one year for the first 104 iterations, then to two years for the next
# 104 iterations, and so on. We stop after 200 iterations. A real run would use
# many more iterations, or a different stopping rule; see
# [How to choose a stopping rule](@ref).

JADE.train(model; iteration_limit = 200, log_frequency = 10)

# ## Simulating the historical years

# `JADE.simulate` wraps `SDDP.simulate` to record the reservoir storage, the
# thermal generation, the price of energy at each node, and the marginal value
# of water in each reservoir.
#
# The default sampling scheme of `SDDP.simulate` samples the inflows of each
# week independently. To see how the policy performs on the historical inflow
# sequences instead, JADE provides two sampling schemes. The first,
# `JADE.HistoricalYears`, has one scenario of 52 weeks for each historical
# year. The scenarios are sampled in order, so we simulate one replication for
# each:

sampling_scheme = JADE.HistoricalYears(data, config)
simulations =
    JADE.simulate(model, length(sampling_scheme.scenarios); sampling_scheme);
length(simulations)

# Each simulation is a vector with one dictionary for each week. The storage is
# a state variable, so it has an incoming and an outgoing value:

simulations[1][1][:x_storage]

# Let's plot the total storage over the year. Each line is one historical
# year. The three years that the policy was trained on are highlighted; the
# other 45 are out-of-sample:

function total_storage(simulation)
    return [sum(x.out for x in week[:x_storage]) for week in simulation]
end

in_sample = [year in config.inflow_years for year in years]
plt = Plots.plot(;
    xlabel = "Week",
    ylabel = "Total storage [Mm³]",
    legend = :bottomright,
)
for (i, simulation) in enumerate(simulations)
    Plots.plot!(
        plt,
        total_storage(simulation);
        color = in_sample[i] ? :red : :gray,
        alpha = in_sample[i] ? 1.0 : 0.4,
        label = i == findfirst(in_sample) ? "In-sample years" :
                i == findfirst(.!in_sample) ? "Out-of-sample years" : "",
    )
end
plt

# The prices and the marginal water values are recorded in \$/MWh. A
# [`SDDP.publication_plot`](@ref) shows the distribution across the 48 years
# of the price at the Haywards node in the peak load block, and of the
# marginal value of water in Lake Pukaki:

Plots.plot(
    SDDP.publication_plot(
        simulations;
        title = "Price at HAY (peak)",
        ylabel = "\$/MWh",
    ) do week
        return week[:price]["HAY", "PEAK"]
    end,
    SDDP.publication_plot(
        simulations;
        title = "Water value of Lake Pukaki",
        ylabel = "\$/MWh",
    ) do week
        return week[:marginal_water_value]["Lake_Pukaki"]
    end;
    xlabel = "Week",
    layout = (1, 2),
    size = (800, 300),
)

# ## Simulating the historical sequence

# The second sampling scheme, `JADE.HistoricalSequence`, has a single scenario
# that runs through the whole historical record, from week 1 of 1970 to week
# 52 of 2017. It shows how the reservoirs would have been managed over the 48
# years if the policy had been in place:

sampling_scheme = JADE.HistoricalSequence(data, config)
simulation = only(JADE.simulate(model, 1; sampling_scheme));
length(simulation)

# The storage, the price, and the water value over the 48 years are:

t = first(years) .+ (0:(length(simulation)-1)) ./ 52
Plots.plot(
    Plots.plot(t, total_storage(simulation); ylabel = "Storage [Mm³]"),
    Plots.plot(
        t,
        [week[:price]["HAY", "PEAK"] for week in simulation];
        ylabel = "Price [\$/MWh]",
    ),
    Plots.plot(
        t,
        [week[:marginal_water_value]["Lake_Pukaki"] for week in simulation];
        ylabel = "Water value [\$/MWh]",
    );
    xlabel = "Year",
    layout = (3, 1),
    legend = false,
    size = (800, 600),
)

# The policy was trained on only three years of inflows and 200 iterations, so
# it should not be read as a good policy for the New Zealand system. The
# spikes in the price are weeks in which the policy sheds load. To improve
# the policy, train on every year and for many more iterations.
