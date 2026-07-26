
# Copyright (c) 2017-26: Oscar Dowson and SDDP.jl contributors.          #src
#                                                                        #src
# This Source Code Form is subject to the terms of the Mozilla Public    #src
# License, v2.0. If a copy of the MPL was not distributed with this file #src
# You can obtain one at http://mozilla.org/MPL/2.0/.                     #src

# # Example: capacity expansion models

# The purpose of this tutorial is to demonstrate a variety of capacity expansion
# problems. The models are intentionally provided with minimal description, you
# should be able to work out what is happening by studying the code.

# This tutorial is an extension of the [Example: deterministic to stochastic](@ref)
# tutorial. We recommend you read that first.

# ## Packages

# This tutorial requires the following packages:

using JuMP
using SDDP
import CSV
import DataFrames
import HiGHS
import Plots

# ## Data

# We use the same data from [Example: deterministic to stochastic](@ref).

io = IOBuffer(
    """
    week,inflow,demand,cost
    1,3,7,10.2\n2,2,7.1,10.4\n3,3,7.2,10.6\n4,2,7.3,10.9\n5,3,7.4,11.2\n
    6,2,7.6,11.5\n7,3,7.8,11.9\n8,2,8.1,12.3\n9,3,8.3,12.7\n10,2,8.6,13.1\n
    11,3,8.9,13.6\n12,2,9.2,14\n13,3,9.5,14.5\n14,2,9.8,14.9\n15,3,10.1,15.3\n
    16,2,10.4,15.8\n17,3,10.7,16.2\n18,2,10.9,16.6\n19,3,11.2,17\n20,3,11.4,17.4\n
    21,3,11.6,17.7\n22,2,11.7,18\n23,3,11.8,18.3\n24,2,11.9,18.5\n25,3,12,18.7\n
    26,2,12,18.9\n27,3,12,19\n28,2,11.9,19.1\n29,3,11.8,19.2\n30,2,11.7,19.2\n
    31,3,11.6,19.2\n32,2,11.4,19.2\n33,3,11.2,19.1\n34,2,10.9,19\n35,3,10.7,18.9\n
    36,2,10.4,18.8\n37,3,10.1,18.6\n38,2,9.8,18.5\n39,3,9.5,18.4\n40,3,9.2,18.2\n
    41,2,8.9,18.1\n42,3,8.6,17.9\n43,2,8.3,17.8\n44,3,8.1,17.7\n45,2,7.8,17.6\n
    46,3,7.6,17.5\n47,2,7.4,17.5\n48,3,7.3,17.5\n49,2,7.2,17.5\n50,3,7.1,17.6\n
    51,3,7,17.7\n52,3,7,17.8\n
    """,
)
data = CSV.read(io, DataFrames.DataFrame)
T = size(data, 1)
reservoir_max = 350.0
reservoir_initial = 300
flow_max = 9

# ## The operational model

# We use the operational problem from [Example: deterministic to stochastic](@ref).

model = SDDP.LinearPolicyGraph(;
    stages = T,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(
        sp,
        0 <= x_storage <= reservoir_max,
        SDDP.State,
        initial_value = reservoir_initial,
    )
    @variable(sp, 0 <= u_flow <= flow_max)
    @variable(sp, 0 <= u_thermal)
    @variable(sp, 0 <= u_spill)
    @variable(sp, ω_inflow)
    Ω, P = [-2, 0, 5], [0.3, 0.4, 0.3]
    SDDP.parameterize(sp, Ω, P) do ω
        fix(ω_inflow, data[t, :inflow] + ω)
        return
    end
    @constraint(sp, x_storage.out == x_storage.in - u_flow - u_spill + ω_inflow)
    @constraint(sp, u_flow + u_thermal == data[t, :demand])
    @stageobjective(sp, data[t, :cost] * u_thermal)
    return
end
SDDP.train(model; iteration_limit = 100)
simulations = SDDP.simulate(model, 100, [:x_storage, :u_flow])
Plots.plot(
    SDDP.publication_plot(simulations; ylabel = "Storage") do sim
        return sim[:x_storage].out
    end,
    SDDP.publication_plot(simulations; ylabel = "Hydro") do sim
        return sim[:u_flow]
    end;
    layout = (2, 1),
)

# ## Invest then operate

# Our operational model assumes a fixed `reservoir_max`. Let's change this to an
# investment decision that we make before we operate the system.

model = SDDP.LinearPolicyGraph(;
    stages = T + 1,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, node
    @variable(sp, x_reservoir_max >= 0, SDDP.State, initial_value = 0)
    @variable(sp, x_storage >= 0, SDDP.State, initial_value = 0)
    @constraint(sp, x_storage.out <= x_reservoir_max.out)
    @variable(sp, 0 <= u_flow <= flow_max)
    @variable(sp, 0 <= u_thermal)
    @variable(sp, 0 <= u_spill)
    @variable(sp, ω_inflow)
    if node == 1  # Investment node
        @stageobjective(sp, x_reservoir_max.out)
        @constraint(sp, x_storage.out <= reservoir_initial)
    else  # Operational node
        t = mod(node - 1, T + 1)
        @constraint(sp, x_reservoir_max.out == x_reservoir_max.in)
        Ω, P = [-2, 0, 5], [0.3, 0.4, 0.3]
        SDDP.parameterize(sp, Ω, P) do ω
            fix(ω_inflow, data[t, :inflow] + ω)
            return
        end
        @constraint(
            sp,
            x_storage.out == x_storage.in - u_flow - u_spill + ω_inflow
        )
        @constraint(sp, u_flow + u_thermal == data[t, :demand])
        @stageobjective(sp, data[t, :cost] * u_thermal)
    end
    return
end
SDDP.train(model; iteration_limit = 100)
simulations = SDDP.simulate(model, 100, [:x_storage, :u_flow, :x_reservoir_max])
Plots.plot(
    SDDP.publication_plot(simulations; ylabel = "Storage") do sim
        return sim[:x_storage].out
    end,
    SDDP.publication_plot(simulations; ylabel = "Hydro") do sim
        return sim[:u_flow]
    end,
    SDDP.publication_plot(simulations; ylabel = "Reservoir Max") do sim
        return sim[:x_reservoir_max].out
    end;
    layout = (2, 2),
)

# ## Multiple investments

# We can also make the `flow_max` an investment option. In general, any bound
# can be made into an investment decision. Add a new state variable. Penalize
# the cost of investment in a new node, and write appropriate constraints for
# the dynamics of the other state variables.

model = SDDP.LinearPolicyGraph(;
    stages = T + 1,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, node
    @variable(sp, x_reservoir_max >= 0, SDDP.State, initial_value = 0)
    @variable(sp, x_flow_max >= 0, SDDP.State, initial_value = 0)
    @variable(sp, x_storage >= 0, SDDP.State, initial_value = 0)
    @constraint(sp, x_storage.out <= x_reservoir_max.out)
    @variable(sp, 0 <= u_flow)
    @constraint(sp, u_flow <= x_flow_max.out)
    @variable(sp, 0 <= u_thermal)
    @variable(sp, 0 <= u_spill)
    @variable(sp, ω_inflow)
    if node == 1  # Investment node
        @stageobjective(sp, x_reservoir_max.out + x_flow_max.out)
        @constraint(sp, x_storage.out <= reservoir_initial)
    else  # Operational node
        t = mod(node - 1, T + 1)
        @constraint(sp, x_reservoir_max.out == x_reservoir_max.in)
        @constraint(sp, x_flow_max.out == x_flow_max.in)
        Ω, P = [-2, 0, 5], [0.3, 0.4, 0.3]
        SDDP.parameterize(sp, Ω, P) do ω
            fix(ω_inflow, data[t, :inflow] + ω)
            return
        end
        @constraint(
            sp,
            x_storage.out == x_storage.in - u_flow - u_spill + ω_inflow
        )
        @constraint(sp, u_flow + u_thermal == data[t, :demand])
        @stageobjective(sp, data[t, :cost] * u_thermal)
    end
    return
end
SDDP.train(model; iteration_limit = 100)
simulations = SDDP.simulate(
    model,
    100,
    [:x_storage, :u_flow, :x_reservoir_max, :x_flow_max],
)
Plots.plot(
    SDDP.publication_plot(simulations; ylabel = "Storage") do sim
        return sim[:x_storage].out
    end,
    SDDP.publication_plot(simulations; ylabel = "Hydro") do sim
        return sim[:u_flow]
    end,
    SDDP.publication_plot(simulations; ylabel = "Reservoir Max") do sim
        return sim[:x_reservoir_max].out
    end,
    SDDP.publication_plot(simulations; ylabel = "Flow Max") do sim
        return sim[:x_flow_max].out
    end;
    layout = (2, 2),
)

# ## Invest-operate-invest-operate

# Now we present a model where we get to change our investments after one year.
# We now have two operational years, each of which is preceeded by an investment
# node. We need to think carefuly about the dynamics and cost of the second
# investment node.

model = SDDP.LinearPolicyGraph(;
    stages = 2 * T + 2,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, node
    @variable(sp, x_reservoir_max >= 0, SDDP.State, initial_value = 0)
    @variable(sp, x_flow_max >= 0, SDDP.State, initial_value = 0)
    @variable(sp, x_storage >= 0, SDDP.State, initial_value = 0)
    @constraint(sp, x_storage.out <= x_reservoir_max.out)
    @variable(sp, 0 <= u_flow)
    @constraint(sp, u_flow <= x_flow_max.out)
    @variable(sp, 0 <= u_thermal)
    @variable(sp, 0 <= u_spill)
    @variable(sp, ω_inflow)
    if node == 1  # First investment node
        @stageobjective(sp, x_reservoir_max.out + x_flow_max.out)
        @constraint(sp, x_storage.out <= reservoir_initial)
    elseif node == T + 2  # Second investment node
        @stageobjective(
            sp,
            (x_reservoir_max.out - x_reservoir_max.in) +
            (x_flow_max.out - x_flow_max.in),
        )
        @constraint(sp, x_storage.out == x_storage.in)
    else  # Operational node
        t = mod(node - 1, T + 1)
        @constraint(sp, x_reservoir_max.out == x_reservoir_max.in)
        @constraint(sp, x_flow_max.out == x_flow_max.in)
        Ω, P = [-2, 0, 5], [0.3, 0.4, 0.3]
        SDDP.parameterize(sp, Ω, P) do ω
            fix(ω_inflow, data[t, :inflow] + ω)
            return
        end
        @constraint(
            sp,
            x_storage.out == x_storage.in - u_flow - u_spill + ω_inflow
        )
        @constraint(sp, u_flow + u_thermal == data[t, :demand])
        @stageobjective(sp, data[t, :cost] * u_thermal)
    end
    return
end
SDDP.train(model; iteration_limit = 100)
simulations = SDDP.simulate(
    model,
    100,
    [:x_storage, :u_flow, :x_reservoir_max, :x_flow_max],
)
Plots.plot(
    SDDP.publication_plot(simulations; ylabel = "Storage") do sim
        return sim[:x_storage].out
    end,
    SDDP.publication_plot(simulations; ylabel = "Hydro") do sim
        return sim[:u_flow]
    end,
    SDDP.publication_plot(simulations; ylabel = "Reservoir Max") do sim
        return sim[:x_reservoir_max].out
    end,
    SDDP.publication_plot(simulations; ylabel = "Flow Max") do sim
        return sim[:x_flow_max].out
    end;
    layout = (2, 2),
)

# ## Invest-operate-invest-operate-loop

# Now we present a model where we enter an infinite horizion operational problem
# after being able to update our investments the second time. For this, we need
# a specialized policy graph:

graph = SDDP.LinearGraph(2 * T + 2)
SDDP.add_edge(graph, 2 * T + 2 => T + 3, 0.95)
model = SDDP.PolicyGraph(
    graph;
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, node
    @variable(sp, x_reservoir_max >= 0, SDDP.State, initial_value = 0)
    @variable(sp, x_flow_max >= 0, SDDP.State, initial_value = 0)
    @variable(sp, x_storage >= 0, SDDP.State, initial_value = 0)
    @constraint(sp, x_storage.out <= x_reservoir_max.out)
    @variable(sp, 0 <= u_flow)
    @constraint(sp, u_flow <= x_flow_max.out)
    @variable(sp, 0 <= u_thermal)
    @variable(sp, 0 <= u_spill)
    @variable(sp, ω_inflow)
    if node == 1  # First investment node
        @stageobjective(sp, x_reservoir_max.out + x_flow_max.out)
        @constraint(sp, x_storage.out <= reservoir_initial)
    elseif node == T + 2  # Second investment node
        @stageobjective(
            sp,
            (x_reservoir_max.out - x_reservoir_max.in) +
            (x_flow_max.out - x_flow_max.in),
        )
        @constraint(sp, x_storage.out == x_storage.in)
    else  # Operational node
        t = mod(node - 1, T + 1)
        @constraint(sp, x_reservoir_max.out == x_reservoir_max.in)
        @constraint(sp, x_flow_max.out == x_flow_max.in)
        Ω, P = [-2, 0, 5], [0.3, 0.4, 0.3]
        SDDP.parameterize(sp, Ω, P) do ω
            fix(ω_inflow, data[t, :inflow] + ω)
            return
        end
        @constraint(
            sp,
            x_storage.out == x_storage.in - u_flow - u_spill + ω_inflow
        )
        @constraint(sp, u_flow + u_thermal == data[t, :demand])
        @stageobjective(sp, data[t, :cost] * u_thermal)
    end
    return
end
SDDP.train(model; iteration_limit = 100)
simulations = SDDP.simulate(
    model,
    100,
    [:x_storage, :u_flow, :x_reservoir_max, :x_flow_max];
    sampling_scheme = SDDP.InSampleMonteCarlo(;
        max_depth = 5 * T + 2,
        terminate_on_dummy_leaf = false,
    ),
)
Plots.plot(
    SDDP.publication_plot(simulations; ylabel = "Storage") do sim
        return sim[:x_storage].out
    end,
    SDDP.publication_plot(simulations; ylabel = "Hydro") do sim
        return sim[:u_flow]
    end,
    SDDP.publication_plot(simulations; ylabel = "Reservoir Max") do sim
        return sim[:x_reservoir_max].out
    end,
    SDDP.publication_plot(simulations; ylabel = "Flow Max") do sim
        return sim[:x_flow_max].out
    end;
    layout = (2, 2),
)
