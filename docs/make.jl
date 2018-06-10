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
            "tutorial/08_odds_and_ends.md",
            "tutorial/09_nonlinear.md",
            "tutorial/10_parallel.md",
            "tutorial/12_price_interpolation.md"
        ],
        "Readings" => "readings.md",
        "Reference" => "apireference.md"
    ],
    assets = [
        "assets/clicked_trajectory.png",
        "assets/dowson_thesis.pdf",
        "assets/plot_value_function.png",
        "assets/publication_plot.png",
        "assets/saddle_function.gif",
        "assets/single_trajectory.png",
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
