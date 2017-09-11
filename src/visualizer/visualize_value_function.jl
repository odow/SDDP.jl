#=
    given a set of cuts, we can fix all but two, discretise those, and evaluate
    the value function the discretized points. Then we can plot the (rough)
    value function. A better way would be to find the facets of the value
    function but I can't get JuliaPolhedra/Polyhedra.jl to install and work
    due to some dependency issues on Windows.
=#

const PLOTLY_HTML_FILE = "valuefunction.template.html"
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

"""
     SDDP.plotvaluefunction(m::SDDPModel, stage::Int, markovstate::Int, states::Union{Float64, AbstractVector{Float64}}...; label1="State 1", label2="State 2")

# Description

Plot the value function of stage `stage` and Markov state `markovstate` in the
SDDPModel `m` at the points in the discretized state space given by `states`. If
the value in `states` is a real number, the state is evaluated at that point. If
the value is a vector, the state is evaluated at all the points in the vector.
At most two states can be vectors.

# Examples

    SDDP.plotvaluefunction(m, 2, 1, 0.0:0.1:1.0, 0.5, 0.0:0.1:1.0; label1="State 1", label2="State 3")
"""
function plotvaluefunction(m::SDDPModel, stage::Int, markovstate::Int, states::Union{Float64, AbstractVector{Float64}}...; label1="State 1", label2="State 2")
    sp = SDDP.getsubproblem(m, stage, markovstate)
    vf = SDDP.valueoracle(sp)
    plotvaluefunction(vf, states...; label1=label1, label2=label2)
end

"""
     SDDP.plotvaluefunction(vf::DefaultValueFunction, states::Union{Float64, AbstractVector{Float64}}...; label1="State 1", label2="State 2")

# Description

Plot the value function `vf` at the points in the discretized state space given
by `states`. If the value in `states` is a real number, the state is evaluated
at that point. If the value is a vector, the state is evaluated at all the points
in the vector. At most two states can be vectors.

# Examples

    SDDP.plotvaluefunction(vf, 0.0:0.1:1.0, 0.5, 0.0:0.1:1.0; label1="State 1", label2="State 3")
"""
function plotvaluefunction(vf::DefaultValueFunction, states::Union{Float64, AbstractVector{Float64}}...; label1="State 1", label2="State 2")
    cuts = SDDP.validcuts(SDDP.cutoracle(vf))
    A, b = getAb(cuts)
    @assert length(states) == size(A, 2)
    free_args = Int[]
    for (i, state) in enumerate(states)
        if isa(state, Float64)
            b += A[:, i] * state
        else
            push!(free_args, i)
        end
    end
    @assert length(free_args) == 2
    A = A[:, free_args]
    x = mesh(states[free_args[1]], states[free_args[2]])
    y = b * ones(1, size(x, 2)) + A * x
    yi = Float64[maximum(y[:, i]) for i in 1:size(y, 2)]::Vector{Float64}
    xdata = x[1,:]
    ydata = x[2,:]
    zdata = yi
    html_string = gethtmlstring(PLOTLY_HTML_FILE)
    html_string = replace(html_string, "<!--xdata-->", json(xdata))
    html_string = replace(html_string, "<!--ydata-->", json(ydata))
    html_string = replace(html_string, "<!--zdata-->", json(zdata))
    html_string = replace(html_string, "<!--xtitle-->", label1)
    html_string = replace(html_string, "<!--ytitle-->", label2)
    launch_file(html_string, PLOTLY_ASSETS)
end
