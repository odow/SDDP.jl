#=
    given a set of cuts, we can fix all but two, discretise those, and evaluate
    the value function the discretized points. Then we can plot the (rough)
    value function. A better way would be to find the facets of the value
    function but I can't get JuliaPolhedra/Polyhedra.jl to install and work
    due to some dependency issues on Windows.
=#
const PLOTLY_ASSET_DIR = dirname(@__FILE__)
const PLOTLY_HTML_FILE = joinpath(PLOTLY_ASSET_DIR, "plotly.valuefunction.template.html")
const PLOTLY_ASSETS    = ["plotly.v1.30.0.js"]

function getAb(cuts::Vector{Cut})
    b = zeros(length(cuts))
    A = zeros(length(cuts), length(cuts[1].coefficients))
    for (i, cut) in enumerate(cuts)
        b[i] = cut.intercept
        for (j, coef) in enumerate(cut.coefficients)
            A[i, j] = coef
        end
    end
    return A, b
end

function mesh(x, y)
    A = zeros(2, length(x) * length(y))
    i = 1
    for xi in x
        for yi in y
            A[1, i] = xi
            A[2, i] = yi
            i += 1
        end
    end
    A
end

function visualizevaluefunction(m::SDDPModel, stage::Int, markovstate::Int, args::Union{Float64, AbstractVector{Float64}}...; state1="State 1", state2="State 2")
    sp = SDDP.getsubproblem(m, stage, markovstate)
    cuts = SDDP.validcuts(SDDP.cutoracle(sp))
    A, b = getAb(cuts)
    @assert length(args) == size(A, 2)
    free_args = Int[]
    for (i, arg) in enumerate(args)
        if isa(arg, Float64)
            b += A[:, i] * arg
        else
            push!(free_args, i)
        end
    end
    @assert length(free_args) == 2
    A = A[:, free_args]
    x = mesh(args[free_args[1]], args[free_args[2]])
    y = b * ones(1, size(x, 2)) + A * x
    yi = Float64[maximum(y[:, i]) for i in 1:size(y, 2)]::Vector{Float64}
    xdata = x[1,:]
    ydata = x[2,:]
    zdata = yi
    html_string = readstring(PLOTLY_HTML_FILE)
    html_string = replace(html_string, "<!--xdata-->", json(xdata))
    html_string = replace(html_string, "<!--ydata-->", json(ydata))
    html_string = replace(html_string, "<!--zdata-->", json(zdata))
    html_string = replace(html_string, "<!--xtitle-->", state1)
    html_string = replace(html_string, "<!--ytitle-->", state2)
    launch_file(html_string, PLOTLY_ASSETS, PLOT_ASSET_DIR)
end
