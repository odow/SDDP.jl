using Documenter, Kokako

makedocs(
    sitename = "Kokako.jl",
    authors  = "Oscar Dowson",
    clean = true,
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    strict = true,
    pages = [
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorial/01_first_steps.md",
            "tutorial/02_adding_uncertainty.md",
            "tutorial/03_objective_uncertainty.md",
            "tutorial/04_markov_uncertainty.md",
            "tutorial/05_plotting.md",
            "tutorial/06_warnings.md",
            "tutorial/11_risk.md",
            "tutorial/12_stopping_rules.md",
            "tutorial/13_generic_graphs.md",
            "tutorial/14_objective_states.md",
            "tutorial/15_performance.md"
        ],
        "Reference" => "apireference.md"
    ],
    assets = [
        "deterministic_linear_policy_graph.png",
        "logo.ico",
        "logo.png",
        "publication_plot.png",
        "spaghetti_plot.html",
        "stochastic_linear_policy_graph.png",
        "stochastic_markovian_policy_graph.png"
    ]
)

deploydocs(repo   = "github.com/odow/Kokako.jl.git")
