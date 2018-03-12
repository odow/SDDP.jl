#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using RecipesBase

function preparedata(s::Vector{Dict{Symbol, Any}}, f::Function, quantiles::Vector{Float64})
    T = length(f(first(s)))
    y = zeros(T, length(s))
    for (i, sd) in enumerate(s)
        y[:,i] .= f(sd)
    end
    hcat([quantile(y[i,:], quantiles) for i in 1:size(y, 1)]...)
end
preparedata(s::Vector{Dict{Symbol, Any}}, key::Symbol, quantiles::Vector{Float64}) = preparedata(s, x->x[key], quantiles)

@userplot PublicationPlot

@recipe function f(h::PublicationPlot; quantile=[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
    dataset, func = h.args
    xlabel --> "Stage\n"
    size --> (500, 300)
    Q = preparedata(dataset, func, sort(quantile))
    for i in 1:floor(Int, size(Q, 1) / 2)
        @series begin
            x := 1:size(Q, 2)
            μ = 0.5 * (Q[i,:] + Q[end- i + 1, :])
            r = Q[end- i + 1, :] - μ
            ribbon := r
            y := μ
            fillalpha--> 0.2
            alpha --> 0.0
            c --> "#00467F"
            label := ""
            ()
        end
    end
    if mod(size(Q, 1) , 2) == 1
        qi = ceil(Int, size(Q, 1) / 2)
        @series begin
            x := 1:size(Q, 2)
            y := Q[qi, :]
            c --> "#00467F"
            label := ""
            ()
        end
    end
end
