var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SDDP-1",
    "page": "Introduction",
    "title": "SDDP",
    "category": "section",
    "text": "(Image: Build Status)"
},

{
    "location": "index.html#Installation-1",
    "page": "Introduction",
    "title": "Installation",
    "category": "section",
    "text": "This package is unregistered so you will need to Pkg.clone it as follows:Pkg.clone(\"https://github.com/odow/SDDP.jl.git\")"
},

{
    "location": "index.html#Note-1",
    "page": "Introduction",
    "title": "Note",
    "category": "section",
    "text": "The documentation is still very incomplete, and the internals of the library need a tidy and a refactor, however the user-facing API from the examples should be stable enough to use."
},

{
    "location": "quick.html#",
    "page": "Quick Start",
    "title": "Quick Start",
    "category": "page",
    "text": ""
},

{
    "location": "quick.html#Quick-Start-1",
    "page": "Quick Start",
    "title": "Quick Start",
    "category": "section",
    "text": "CurrentModule = SDDP"
},

{
    "location": "quick.html#SDDP.SDDPModel",
    "page": "Quick Start",
    "title": "SDDP.SDDPModel",
    "category": "Type",
    "text": "SDDPModel(;kwargs...) do ...\n\nend\n\nDescription\n\nThis function constructs an SDDPModel.\n\nRequired Keyword arguments\n\nstages::Int\n\nThe number of stages in the problem. A stage is defined as each step in time at  which a decion can be made. Defaults to 1.\n\nobjective_bound::Float64\nsolver::MathProgBase.AbstractMathProgSolver\n\nOptional Keyword arguments\n\ncut_oracle\nrisk_measure\nnoise_probability\nmarkov_transition\n\nReturns\n\n* `m`: the `SDDPModel`\n\n\n\n"
},

{
    "location": "quick.html#Initialising-the-model-object-1",
    "page": "Quick Start",
    "title": "Initialising the model object",
    "category": "section",
    "text": "The first step is to initialise the SDDP model object. We do this using the following syntax:If we have more than one markov state:m = SDDPModel([;kwargs...]) do sp, stage, markov_state\n\n    # Stage problem definition where `sp` is a `JuMP.Model`object,\n\nendOtherwise if we have a single markov statem = SDDPModel([;kwargs...]) do sp, stage\n\n  # Stage problem definition\n\nendSDDPModel"
},

{
    "location": "quick.html#JuMP.solve",
    "page": "Quick Start",
    "title": "JuMP.solve",
    "category": "Function",
    "text": "solve(m::SDDPModel; kwargs...)\n\nDescription\n\nSolve the SDDPModel m using SDDP. Accepts a number of keyword arguments to control the solution process.\n\nPositional arguments\n\nm: the SDDPModel to solve\n\nKeyword arguments\n\nmax_iterations::Int:  The maximum number of cuts to add to a single stage problem before terminating.  Defaults to 10.\ntime_limit::Real:  The maximum number of seconds (in real time) to compute for before termination.  Defaults to Inf.\nsimulation::MonteCarloSimulation:  We control the behaviour of the policy simulation phase of the algorithm using  the MonteCarloSimulation(;kwargs...) constructor. This just groups a  series of related keyword arguments. The keywords are\nfrequency::Int\nThe frequency (by iteration) with which to run the policy simulation phase of  the algorithm in order to construct a statistical bound for the policy. Defaults  to 0 (never run).\nmin::Float64\nMinimum number of simulations to conduct before constructing a confidence interval  for the bound. Defaults to 20.\nstep::Float64\nNumber of additional simulations to conduct before constructing a new confidence  interval for the bound. Defaults to 1.\nmax::Float64\nMaximum number of simulations to conduct in the policy simulation phase. Defaults  to min.\nconfidence::Float64\nConfidence level of the confidence interval. Defaults to 0.95 (95% CI).\ntermination::Bool\nWhether to terminate the solution algorithm with the status :converged if the  deterministic bound is with in the statistical bound after max simulations.  Defaults to false.\nbound_convergence:  We may also wish to terminate the algorithm if the deterministic bound stalls  for a specified number of iterations (regardless of whether the policy has  converged). This can be controlled by the BoundConvergence(;kwargs...)  constructor. It has the following keywords:\niterations::Int\nTerminate if the maximum deviation in the deterministic bound from the mean  over the last iterations number of iterations is less than rtol (in  relative terms) or atol (in absolute terms).\nrtol::Float64\nMaximum allowed relative deviation from the mean.  Defaults to 0.0\natol::Float64\nMaximum allowed absolute deviation from the mean.  Defaults to 0.0\ncut_selection_frequency::Int:  Frequency (by iteration) with which to rebuild subproblems using a subset of  cuts. Frequent cut selection (i.e. cut_selection_frequency is small) reduces  the size of the subproblems that are solved, but incurrs the overhead of rebuilding  the subproblems. However, infrequent cut selection (i.e.  cut_selection_frequency is large) allows the subproblems to grow large (many  constraints) leading to an increase in the solve time of individual subproblems.  Defaults to 0 (never run).\nprint_level::Int:   0 - off: nothing logged to screen (still written to log file if specified).   1 - on: solve iterations written to screen.   Defaults to 1\nlog_file::String:  Relative filename to write the log to disk. Defaults to \"\" (no log written)\nsolve_type:  One of\nAsyncronous() - solve using a parallelised algorithm\nSerial() - solve using a serial algorithm\nDefault chooses automatically based on the number of available processors.\nreduce_memory_footprint::Bool:  Implements the idea proposed in https://github.com/JuliaOpt/JuMP.jl/issues/969#issuecomment-282191105  to reduce the memory consumption when running SDDP. This is an issue if you  wish to save the model m to disk since it discards important information.  Defaults to false.\ncut_output_file::String:  Relative filename to write discovered cuts to disk. Defaults to \"\" (no cuts written)\n\nReturns\n\nstatus::Symbol:  Reason for termination. One of\n:solving\n:interrupted\n:converged\n:max_iterations\n:bound_convergence\n:time_limit\n\n\n\n"
},

{
    "location": "quick.html#Solve-1",
    "page": "Quick Start",
    "title": "Solve",
    "category": "section",
    "text": "To solve the SDDP model m we use status = solve(m::SDDPModel; kwargs...). This accepts a number of keyword arguments to control the solution process.solve"
},

{
    "location": "quick.html#Simulate-1",
    "page": "Quick Start",
    "title": "Simulate",
    "category": "section",
    "text": ""
},

{
    "location": "quick.html#SDDP.@visualise",
    "page": "Quick Start",
    "title": "SDDP.@visualise",
    "category": "Macro",
    "text": "@visualise(results, replication, stage, begin\n	... plot definitions ...\nend)\n\nDescription\n\nPlot everything using interactive javascript. This will launch an HTML page\nto explore.\n\nUsage\n\n    `@visualise(results, i, t, begin\n        ... one line for each plot ...\n    end)`\n\nwhere results is the vector of result dictionaries from simulate(), i is the\nsimulation index (1:length(results)), and t is the stage index (1:T).\n\nEach plot line gets transformed into an anonymous function\n    (results, i, t) -> ... plot line ...\nso can be any valid Julia syntax that uses results, i, or t as an argument.\n\nAfter the plot definition, keyword arguments can be used (in parenthesises):\n    `title`       - set the title of the plot\n    `ylabel`      - set the yaxis label\n    `xlabel`      - set the xaxis label\n    `interpolate` - interpolate lines between stages. Defaults to \"linear\"\n    see https://github.com/d3/d3-3.x-api-reference/blob/master/SVG-Shapes.md\n        #line_interpolate for all options\n\nResults Object\n\n`results::Vector{Dict{Symbol, Any}}` is a vector of dictionaries where each\ndictionary corresponds to one simulation (therefore there will be\n`N = length(results)` lines plotted in each graph).\n\n\n\n"
},

{
    "location": "quick.html#Visualise-1",
    "page": "Quick Start",
    "title": "Visualise",
    "category": "section",
    "text": "@visualise"
},

]}
