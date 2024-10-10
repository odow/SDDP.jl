#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Internal function: convert dataset (from SDDP.simulate) into a matrix where
# the rows are quantiles, and the columns are stages.
function publication_data(
    dataset::Vector{<:Vector{<:AbstractDict}},
    quantiles::Vector{Float64},
    stage_function::Function,
)
    max_stages = maximum(length.(dataset))
    output_array = fill(NaN, length(quantiles), max_stages)
    for t in 1:max_stages
        stage_data = Float64[]
        for (i, d) in enumerate(dataset)
            if length(d) < t
                # Skip this replication because it doesn't have enough stages
                continue
            end
            s = stage_function(d[t])
            if !isfinite(s)
                error(
                    "Unable to plot `publication_plot` because stage $t " *
                    "of replication $i contains data that is not finite. " *
                    "The data function must return a finite real-valued " *
                    "scalar. Got: $s",
                )
            end
            push!(stage_data, s)
        end
        output_array[:, t] .= Statistics.quantile(stage_data, quantiles)
    end
    return output_array
end

"""
    SDDP.publication_plot(
        data_function, simulations;
        quantile = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0],
        kwargs...)

Create a `Plots.jl` recipe plot of the simulations.

See `Plots.jl` for the list of keyword arguments.

## Examples

    SDDP.publication_plot(simulations; title = "My title") do data
        return data[:stage_objective]
    end
"""
function publication_plot(
    data_function::Function,
    simulations::Vector{<:Vector{<:AbstractDict}};
    kwargs...,
)
    # An annoying over-load so that we can provide a consistent interface
    # instead of the Plots.jl generated `publicationplot`.
    return SDDP.publicationplot(simulations, data_function; kwargs...)
end

RecipesBase.@userplot PublicationPlot

RecipesBase.@recipe function f(
    publication_plot::PublicationPlot;
    quantile = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0],
)
    dataset, stage_function = publication_plot.args
    size --> (500, 300)
    data_matrix = publication_data(dataset, sort(quantile), stage_function)
    for i in 1:floor(Int, size(data_matrix, 1) / 2)
        μ = 0.5 * (data_matrix[i, :] + data_matrix[end-i+1, :])
        r = data_matrix[end-i+1, :] - μ
        RecipesBase.@series begin
            x := 1:size(data_matrix, 2)
            ribbon := r
            y := μ
            fillalpha --> 0.2
            seriesalpha --> 0.0
            seriescolor --> "#00467F"
            label := ""
            ()
        end
    end
    if mod(size(data_matrix, 1), 2) == 1
        qi = ceil(Int, size(data_matrix, 1) / 2)
        RecipesBase.@series begin
            x := 1:size(data_matrix, 2)
            y := data_matrix[qi, :]
            seriescolor --> "#00467F"
            label := ""
            ()
        end
    end
end
