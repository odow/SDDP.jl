using Documenter, SDDP

makedocs(
    doctest  = false,
    clean    = true,
    format   = :html,
    sitename = "SDDP.jl",
    authors  = "Oscar Dowson",
    pages = [
        "Home" => "index.md",
        # "Introduction" => "introduction.md",
        "Tutorials" => Any[
            "tutorial/01_first_steps.md",
            "tutorial/02_rhs_noise.md",
            "tutorial/03_objective_noise.md",
            "tutorial/04_markovian_policygraphs.md",
            "tutorial/05_risk.md",
            "tutorial/06_cut_selection.md",
            "tutorial/07_plotting.md",
            "tutorial/08_odds_and_ends.md"
        ],
        "Readings" => "readings.md",
        "Old Manual" => "oldindex.md",
        "Reference" => "apireference.md"
    ],
    assets = [
        "clicked_trajectory.png",
        "single_trajectory.png",
        "publication_plot.png",
        "plot_value_function.png"
    ]
)

deploydocs(
    repo   = "github.com/odow/SDDP.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing,
)
