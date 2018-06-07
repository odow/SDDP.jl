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
            "tutorial/06_cut_selection.md"
        ],
        "Readings" => "readings.md",
        "Old Manual" => "oldindex.md",
        "Reference" => "apireference.md"
    ],
    assets = ["3d.gif"]
)

deploydocs(
    repo   = "github.com/odow/SDDP.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing,
)
