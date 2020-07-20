# # Basic V: plotting

# In our previous tutorials, we formulated, solved, and simulated multistage
# stochastic optimization problems. However, we haven't really investigated what
# the solution looks like. Luckily, `SDDP.jl` includes a number of plotting
# tools to help us do that. In this tutorial, we explain the tools and make some
# pretty pictures.

# ## Preliminaries

# The next two plot types help visualize the policy. Thus, we first need to
# create a policy and simulate some trajectories. So, let's take the model from
# [Basic IV: Markov uncertainty](@ref), train it for 20 iterations, and then
# simulate 100 Monte Carlo realizations of the policy.

using SDDP, GLPK

Ω = [
    (inflow = 0.0, fuel_multiplier = 1.5),
    (inflow = 50.0, fuel_multiplier = 1.0),
    (inflow = 100.0, fuel_multiplier = 0.75)
]

model = SDDP.MarkovianPolicyGraph(
    transition_matrices = Array{Float64, 2}[
        [1.0]', [0.75 0.25], [0.75 0.25 ; 0.25 0.75]
    ],
    sense = :Min,
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer
) do subproblem, node
    t, markov_state = node
    @variable(subproblem, 0 <= volume <= 200, SDDP.State, initial_value = 200)
    @variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation >= 0
        hydro_spill >= 0
        inflow
    end)
    @constraints(subproblem, begin
        volume.out == volume.in + inflow - hydro_generation - hydro_spill
        thermal_generation + hydro_generation == 150.0
    end)
    probability = markov_state == 1 ? [1/6, 1/3, 1/2] : [1/2, 1/3, 1/6]
    fuel_cost = [50.0, 100.0, 150.0]
    SDDP.parameterize(subproblem, Ω, probability) do ω
        JuMP.fix(inflow, ω.inflow)
        @stageobjective(subproblem,
            ω.fuel_multiplier * fuel_cost[t] * thermal_generation
        )
    end
end

SDDP.train(model, iteration_limit = 20, run_numerical_stability_report = false)

simulations = SDDP.simulate(
    model,
    100,
    [:volume, :thermal_generation, :hydro_generation, :hydro_spill])

println("Completed $(length(simulations)) simulations.")

# Great! Now we have some data in `simulations` to visualize.

# ## Spaghetti plots

# The first plotting utility we discuss is a _spaghetti_ plot (you'll understand
# the name when you see the graph).

# To create a spaghetti plot, begin by creating a new
# [`SDDP.SpaghettiPlot`](@ref) instance as follows:

plt = SDDP.SpaghettiPlot(simulations)


# We can add plots to `plt` using the [`SDDP.add_spaghetti`](@ref) function.

SDDP.add_spaghetti(plt; title = "Reservoir volume") do data
    return data[:volume].out
end

# You don't have just return values from the simulation, you can compute things
# too.

SDDP.add_spaghetti(
    plt; title = "Fuel cost", ymin = 0, ymax = 250
) do data
    if data[:thermal_generation] > 0
        return data[:stage_objective] / data[:thermal_generation]
    else  # No thermal generation, so return 0.0.
        return 0.0
    end
end

# Note that there are many keyword arguments in addition to `title`. For
# example, we fixed the minimum and maximum values of the y-axis using `ymin`
# and `ymax`. See the [`SDDP.add_spaghetti`](@ref) documentation for all the
# arguments.

# Having built the plot, we now need to display it using [`SDDP.plot`](@ref).

# ```julia
# SDDP.plot(plt, "spaghetti_plot.html")
# ```

# This should open a webpage that looks like [this
# one](../assets/spaghetti_plot.html).

# Using the mouse, you can highlight individual trajectories by hovering over
# them. This makes it possible to visualize a single trajectory across multiple
# dimensions.

# If you click on the plot, then trajectories that are close to the mouse
# pointer are shown darker and those further away are shown lighter.

# ## Publication plots

# Instead of the interactive Javascript plots, you can also create some
# publication ready plots using the [`SDDP.publication_plot`](@ref) function.

# !!! info
#     You need to install the [Plots.jl](https://github.com/JuliaPlots/Plots)
#     package for this to work. We used the `GR` backend (`gr()`), but any
#     `Plots.jl` backend should work.

# [`SDDP.publication_plot`](@ref) implements a plot recipe to create ribbon
# plots of each variable against the stages. The first argument is the vector of
# simulation dictionaries and the second argument is the dictionary key that you
# want to plot. Standard `Plots.jl` keyword arguments such as `title` and `xlabel`
# can be used to modify the look of each plot. By default, the plot displays
# ribbons of the 0-100, 10-90, and 25-75 percentiles. The dark, solid line in the
# middle is the median (i.e. 50'th percentile).

# ```julia
# using Plots
# plot(
#     SDDP.publication_plot(simulations, title = "Outgoing volume") do data
#         return data[:volume].out
#     end,
#     SDDP.publication_plot(simulations, title = "Thermal generation") do data
#         return data[:thermal_generation]
#     end,
#     xlabel = "Stage",
#     ylims = (0, 200),
#     layout = (1, 2),
#     margin_bottom = 5,
#     size = (1000, 300)
# )
# ```

# This should open a plot window with a plot that looks like:
#
# ![publication plot](../assets/publication_plot.png)

# You can save this plot as a PDF using the `Plots.jl` function `savefig`:
# ```julia
# Plots.savefig("my_picture.pdf")
# ```

# ## Plotting the value function

# You can obtain an object representing the value function of a node using
# [`SDDP.ValueFunction`](@ref).

V = SDDP.ValueFunction(model[(1, 1)])

# The value function can be evaluated using [`SDDP.evaluate`](@ref).

SDDP.evaluate(V; volume = 1)

# `evaluate` returns the heigh of the value function, and a subgradient with respect to the
# convex state variables.

# You can also plot the value function using [`SDDP.plot`](@ref)

SDDP.plot(V, volume = 0:200, filename = "value_function.html")

# This should open a webpage that looks like [this one](../assets/value_function.html).

# ## Convergence dashboard

# If the text-based logging isn't to your liking, you can open a visualization of
# the training by passing  `dashboard = true` to [`SDDP.train`](@ref).
# ```julia
# SDDP.train(model; dashboard = true)
# ```
# By default, `dashboard = false` because there is an initial overhead
# associated with opening and preparing the plot.

# !!! warning
#     The dashboard is experimental. There are known bugs associated with it,
#     e.g., [SDDP.jl#226](https://github.com/odow/SDDP.jl/issues/226).

# This concludes our fifth tutorial for `SDDP.jl`. In our next tutorial, [Basic
# VI: words of warning](@ref) we discuss some of the issues that you should be
# aware of when creating your own models.
