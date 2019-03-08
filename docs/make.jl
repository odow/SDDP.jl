using Documenter, Kokako

"Call julia docs/make.jl --fix to rebuild the doctests."
const FIX_DOCTESTS = length(ARGS) == 1 && ARGS[1] == "--fix"

# doctest=:fix only works with `\n` line endings, so let's replace any `\r\n`
# ones.
function convert_line_endings()
    for file in readdir(joinpath(@__DIR__, "src", "tutorial"))
        !endswith(file, ".md") && continue
        filename = joinpath(@__DIR__, "src", "tutorial", file)
        code = read(filename, String)
        open(filename, "w") do io
            write(io, replace(code, "\r\n" => "\n"))
        end
    end
end
FIX_DOCTESTS && convert_line_endings()

makedocs(
    sitename = "SDDP.jl",
    authors  = "Oscar Dowson",
    clean = true,
    doctest = FIX_DOCTESTS ? :fix : true,
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    strict = true,
    pages = [
        "Home" => "index.md",
        "Tutorials" => Any[
            "Basic" => Any[
                "tutorial/01_first_steps.md",
                "tutorial/02_adding_uncertainty.md",
                "tutorial/03_objective_uncertainty.md",
                "tutorial/04_markov_uncertainty.md",
                "tutorial/05_plotting.md",
                "tutorial/06_warnings.md",
                "tutorial/07_advanced_modelling.md"
            ],
            "Intermediate" => Any[
                "tutorial/11_risk.md",
                "tutorial/12_stopping_rules.md",
                "tutorial/13_generic_graphs.md",
                "tutorial/14_objective_states.md",
                "tutorial/15_performance.md"
            ]
        ],
        "Reference" => "apireference.md"
    ],
    assets = [
        "deterministic_linear_policy_graph.png",
        "publication_plot.png",
        "spaghetti_plot.html",
        "stochastic_linear_policy_graph.png",
        "stochastic_markovian_policy_graph.png"
    ],
    doctestfilters = [r"[\s\-]?\d\.\d{6}e[\+\-]\d{2}"]
)

deploydocs(repo = "github.com/odow/Kokako.jl.git")
