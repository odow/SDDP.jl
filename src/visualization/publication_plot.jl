#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Internal function: convert dataset (from SDDP.simulate) into a matrix where
# the rows are quantiles, and the columns are stages.
function publication_data(
    dataset::Vector{Vector{Dict{Symbol,Any}}},
    quantiles::Vector{Float64},
    stage_function::Function,
)
    max_stages = maximum(length.(dataset))
    output_array = fill(NaN, length(quantiles), max_stages)
    for stage = 1:max_stages
        stage_data = stage_function.([data[stage] for data in dataset])
        output_array[:, stage] .= Statistics.quantile(stage_data, quantiles)
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

### Example

    SDDP.publication_plot(simulations; title = "My title") do data
        return data[:stage_objective]
    end
"""
function publication_plot(
    data_function::Function,
    simulations::Vector{Vector{Dict{Symbol,Any}}};
    kwargs...,
)
    # An annoying over-load so that we can provide a consistent interface
    # instead of the Plots.jl generated `publicationplot`.
    return SDDP.publicationplot(simulations, data_function, kwargs...)
end

RecipesBase.@userplot PublicationPlot

RecipesBase.@recipe function f(
    publication_plot::PublicationPlot;
    quantile = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0],
)
    dataset, stage_function = publication_plot.args
    size --> (500, 300)
    data_matrix = publication_data(dataset, sort(quantile), stage_function)
    for i = 1:floor(Int, size(data_matrix, 1) / 2)
        RecipesBase.@series begin
            x := 1:size(data_matrix, 2)
            μ = 0.5 * (data_matrix[i, :] + data_matrix[end-i+1, :])
            r = data_matrix[end-i+1, :] - μ
            ribbon := r
            y := μ
            fillalpha --> 0.2
            alpha --> 0.0
            c --> "#00467F"
            label := ""
            ()
        end
    end
    if mod(size(data_matrix, 1), 2) == 1
        qi = ceil(Int, size(data_matrix, 1) / 2)
        RecipesBase.@series begin
            x := 1:size(data_matrix, 2)
            y := data_matrix[qi, :]
            c --> "#00467F"
            label := ""
            ()
        end
    end
end
