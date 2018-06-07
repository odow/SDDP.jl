using Documenter, SDDP

makedocs(
    doctest  = false,
    clean    = true,
    format   = :html,
    sitename = "SDDP.jl",
    authors  = "Oscar Dowson",
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "Tutorials" => Any[
            "tutorial/first_example.md",
            "tutorial/rhs_noise.md",
            "tutorial/objective_noise.md",
            "tutorial/markovian_policygraphs.md",
            "tutorial/simulation.md",
            "tutorial/risk.md",
            "tutorial/cut_selection.md"
        ],
        "Readings" => "readings.md",
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
