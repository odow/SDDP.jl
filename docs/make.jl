using Documenter, Kokako

makedocs(
    sitename = "SDDP.jl",
    authors  = "Oscar Dowson",
    clean = true,
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/JuliaOpt/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorial/01_first_steps.md",
            "tutorial/02_adding_uncertainty.md",
            "tutorial/03_objective_uncertainty.md",
            "tutorial/04_markov_uncertainty.md"
            # "tutorial/05_risk.md",
            # "tutorial/06_cut_selection.md",
            # "tutorial/07_plotting.md",
            # "tutorial/08_odds_and_ends.md",
            # "tutorial/09_nonlinear.md",
            # "tutorial/10_parallel.md",
            # "tutorial/11_DRO.md",
            # "tutorial/12_price_interpolation.md",
            # "tutorial/13_constraint_noise.md"
        ],
        "Reference" => "apireference.md"
    ],
    assets = [
        "deterministic_linear_policy_graph.png",
        "stochastic_linear_policy_graph.png",
        "stochastic_markovian_policy_graph.png"
    ]
)

deploydocs(repo   = "github.com/odow/Kokako.jl.git")
