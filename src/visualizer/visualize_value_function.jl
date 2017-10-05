#=
    given a set of cuts, we can fix all but two, discretise those, and evaluate
    the value function the discretized points. Then we can plot the (rough)
    value function. A better way would be to find the facets of the value
    function but I can't get JuliaPolhedra/Polyhedra.jl to install and work
    due to some dependency issues on Windows.
=#

const PLOTLY_HTML_FILE = "valuefunction.template.html"
const PLOTLY_ASSETS    = ["plotly.v1.30.0.js"]

const PLOTLY_DATA = Dict{String, Any}(
    # "xaxis" => "x",
    # "yaxis" => "y",
    # "x" => <!--xdata-->,
    # "y" => <!--ydata-->,
    # "z" => <!--zdata-->,
    "showlegend" => false,
    "marker" => Dict(
        "symbol"=>"circle",
        "color"=>"rgba(0, 154, 250, 1.000)",
        "size"=>2
    )
    # "type" => "scatter3d"
)
const PLOTLY_SCENE = Dict{String, Any}()

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

function processvaluefunctiondata(vf::DefaultValueFunction, is_minimization::Bool, states::Union{Float64, AbstractVector{Float64}}...)
    cuts = SDDP.validcuts(SDDP.cutoracle(vf))
    processvaluefunctiondata(cuts, is_minimization, states...)
end
function processvaluefunctiondata(cuts::Vector{Cut}, is_minimization::Bool, states::Union{Float64, AbstractVector{Float64}}...)
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
    @assert length(free_args) == 1 || length(free_args) == 2
    A = A[:, free_args]
    if length(free_args) == 1
        x = states[free_args[1]]'
    else
        x = mesh(states[free_args[1]], states[free_args[2]])
    end
    y = b * ones(1, size(x, 2)) + A * x
    if is_minimization
        yi = Float64[maximum(y[:, i]) for i in 1:size(y, 2)]::Vector{Float64}
    else
        yi = Float64[minimum(y[:, i]) for i in 1:size(y, 2)]::Vector{Float64}
    end
    return x, yi
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
    html = prepvaluefunctionplot(m, stage, markovstate, label1, label2, states...)
    launch_file(html, PLOTLY_ASSETS)
end

function plotvaluefunction(filename::String, m::SDDPModel, stage::Int, markovstate::Int, states::Union{Float64, AbstractVector{Float64}}...; label1="State 1", label2="State 2")
    html = prepvaluefunctionplot(m, stage, markovstate, label1, label2, states...)
    launch_file(html, PLOTLY_ASSETS, filename)
end


function prepvaluefunctionplot(m::SDDPModel, stage::Int, markovstate::Int, label1, label2, states::Union{Float64, AbstractVector{Float64}}...)
    sp = SDDP.getsubproblem(m, stage, markovstate)
    vf = SDDP.valueoracle(sp)
    is_minimization = getsense(sp) == :Min
    x, yi  = processvaluefunctiondata(vf, is_minimization, states...)
    (plotly_data, scene_text) = getplotlydata(x, yi, label1, label2)
    return prephtml(PLOTLY_HTML_FILE, ("<!--DATA-->", json(plotly_data)), ("<!--SCENE-->", scene_text))
end

function getplotlydata(x, yi, label1, label2)
    plotly_data = deepcopy(PLOTLY_DATA)
    plotly_data["x"] = x[1,:]
    if size(x, 1) == 1
        plotly_data["y"] = yi
        plotly_data["type"] = "scatter"
        scene_text = "xaxis: {title: \"$(label1)\"}, yaxis: {title: \"Future Cost\"}"
    else
        plotly_data["xaxis"] = "x"
        plotly_data["yaxis"] = "y"
        plotly_data["type"] = "scatter3d"
        plotly_data["y"] = x[2,:]
        plotly_data["z"] = yi
        plotly_data["mode"] = "markers"
        scene = Dict{String, Any}()
        scene["xaxis"] = Dict("title"=>label1)
        scene["yaxis"] = Dict("title" => label2)
        scene["zaxis"] = Dict("title" => "Future Cost")
        scene["aspectratio"] = Dict("x" => 1, "y" => 1, "z" => 1)
        scene_text = "\"scene\": $(json(scene))"
    end
    return plotly_data, scene_text
end
