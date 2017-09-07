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
    "text": "solve(m::SDDPModel; kwargs...)\n\nDescription\n\nSolve the SDDPModel m using SDDP. Accepts a number of keyword arguments to control the solution process.\n\nPositional arguments\n\nm: the SDDPModel to solve\n\nKeyword arguments\n\nmax_iterations::Int:  The maximum number of cuts to add to a single stage problem before terminating.  Defaults to 10.\ntime_limit::Real:  The maximum number of seconds (in real time) to compute for before termination.  Defaults to Inf.\nsimulation::MonteCarloSimulation: see MonteCarloSimulation\nbound_convergence::BoundConvergence: see BoundConvergence\ncut_selection_frequency::Int:  Frequency (by iteration) with which to rebuild subproblems using a subset of  cuts. Frequent cut selection (i.e. cut_selection_frequency is small) reduces  the size of the subproblems that are solved, but incurrs the overhead of rebuilding  the subproblems. However, infrequent cut selection (i.e.  cut_selection_frequency is large) allows the subproblems to grow large (many  constraints) leading to an increase in the solve time of individual subproblems.  Defaults to 0 (never run).\nprint_level::Int:   0 - off: nothing logged to screen (still written to log file if specified).   1 - on: solve iterations written to screen.   Defaults to 1\nlog_file::String:  Relative filename to write the log to disk. Defaults to \"\" (no log written)\nsolve_type:  One of\nAsyncronous() - solve using a parallelised algorithm\nSerial() - solve using a serial algorithm\nDefault chooses automatically based on the number of available processors.\nreduce_memory_footprint::Bool:  Implements the idea proposed in https://github.com/JuliaOpt/JuMP.jl/issues/969#issuecomment-282191105  to reduce the memory consumption when running SDDP. This is an issue if you  wish to save the model m to disk since it discards important information.  Defaults to false.\ncut_output_file::String:  Relative filename to write discovered cuts to disk. Defaults to \"\" (no cuts written)\n\nReturns\n\nstatus::Symbol:  Reason for termination. One of\n:solving\n:interrupted\n:converged\n:max_iterations\n:bound_convergence\n:time_limit\n\n\n\n"
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

{
    "location": "apireference.html#",
    "page": "Reference",
    "title": "Reference",
    "category": "page",
    "text": "CurrentModule = SDDP"
},

{
    "location": "apireference.html#API-Reference-1",
    "page": "Reference",
    "title": "API Reference",
    "category": "section",
    "text": ""
},

{
    "location": "apireference.html#Defining-the-model-1",
    "page": "Reference",
    "title": "Defining the model",
    "category": "section",
    "text": "SDDPModel"
},

{
    "location": "apireference.html#SDDP.@state",
    "page": "Reference",
    "title": "SDDP.@state",
    "category": "Macro",
    "text": "@state(sp, stateleaving, stateentering)\n\nDescription\n\nDefine a new state variable in the subproblem sp.\n\nArguments\n\nsp               the subproblem\nstateleaving     any valid JuMP @variable syntax to define the value of the state variable at the end of the stage\nstateentering    any valid JuMP @variable syntax to define the value of the state variable at the beginning of the stage\n\nExamples\n\n@state(sp, 0 <= x[i=1:3] <= 1, x0==rand(3)[i] )\n@state(sp,      y        <= 1, y0==0.5        )\n@state(sp,      z            , z0==0.5        )\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.@states",
    "page": "Reference",
    "title": "SDDP.@states",
    "category": "Macro",
    "text": "@states(sp, begin\n    stateleaving1, stateentering1\n    stateleaving2, stateentering2\nend)\n\nDescription\n\nDefine a new state variables in the subproblem sp.\n\nArguments\n\nsp               the subproblem\nstateleaving     any valid JuMP @variable syntax to define the value of the state variable at the end of the stage\nstateentering    any valid JuMP @variable syntax to define the value of the state variable at the beginning of the stage\n\nUsage\n\n@states(sp, begin\n    0 <= x[i=1:3] <= 1, x0==rand(3)[i]\n         y        <= 1, y0==0.5\n         z            , z0==0.5\n end)\n\n\n\n"
},

{
    "location": "apireference.html#States-1",
    "page": "Reference",
    "title": "States",
    "category": "section",
    "text": "@state\n@states"
},

{
    "location": "apireference.html#SDDP.@noise",
    "page": "Reference",
    "title": "SDDP.@noise",
    "category": "Macro",
    "text": "@noise(sp, rhs, constraint)\n\nDescription\n\nAdd a constraint with a noise in the RHS vector to the subproblem sp.\n\nArguments\n\nsp         the subproblem\nrhs        keyword argument key=value where value is a one-dimensional array containing the noise realisations\nconstraint any valid JuMP @constraint syntax that includes the keyword defined by rhs\n\nExamples\n\n@noise(sp, i=1:2, x + y <= i )\n@noise(sp, i=1:2, x + y <= 3 * rand(2)[i] )\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.@noises",
    "page": "Reference",
    "title": "SDDP.@noises",
    "category": "Macro",
    "text": "@noises(sp, rhs, begin\n    constraint\nend)\n\nDescription\n\nThe plural form of @noise similar to the JuMP macro @constraints.\n\nArguments\n\nSee @noise.\n\nExamples\n\n@noises(sp, i=1:2, begin\n    x + y <= i\n    x + y <= 3 * rand(2)[i]\nend)\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.setnoiseprobability!",
    "page": "Reference",
    "title": "SDDP.setnoiseprobability!",
    "category": "Function",
    "text": "setnoiseprobability!(sp::JuMP.Model, distribution::Vector{Float64})\n\nDescription\n\nSet the probability distribution of the stagewise independent noise in the sp subproblem.\n\nArguments\n\nsp            the subproblem\ndistribution vector containing the probability of each outcome occuring.   Should sum to 1. Defaults to the uniform distribution.\n\nExamples\n\nIf there are two realizations:\n\nsetnoiseprobability!(sp, [0.3, 0.7])\nsetnoiseprobability!(sp, [0.5, 0.6]) will error!\n\n\n\n"
},

{
    "location": "apireference.html#Noises-1",
    "page": "Reference",
    "title": "Noises",
    "category": "section",
    "text": "@noise\n@noises\nsetnoiseprobability!"
},

{
    "location": "apireference.html#Objective-1",
    "page": "Reference",
    "title": "Objective",
    "category": "section",
    "text": "stageobjective!"
},

{
    "location": "apireference.html#SDDP.AbstractRiskMeasure",
    "page": "Reference",
    "title": "SDDP.AbstractRiskMeasure",
    "category": "Type",
    "text": "AbstractRiskMeasure\n\nDescription\n\nAbstract type for all risk measures.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.Expectation",
    "page": "Reference",
    "title": "SDDP.Expectation",
    "category": "Type",
    "text": "Expectation()\n\nDescription\n\nThe expectation risk measure.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.NestedAVaR",
    "page": "Reference",
    "title": "SDDP.NestedAVaR",
    "category": "Type",
    "text": "NestedAVaR(;lambda=1.0, beta=1.0)\n\nDescription\n\nA risk measure that is a convex combination of Expectation and Average Value @ Risk (also called Conditional Value @ Risk).\n\nλ * E[x] + (1 - λ) * AV@R(1-β)[x]\n\nKeyword Arguments\n\nlambda\n\nConvex weight on the expectation ((1-lambda) weight is put on the AV@R component. Inreasing values of lambda are less risk averse (more weight on expecattion)\n\nbeta\n\nThe quantile at which to calculate the Average Value @ Risk. Increasing values  of beta are less risk averse. If beta=0, then the AV@R component is the  worst case risk measure.\n\nReturns\n\n`m::NestedAVaR<:AbstractRiskMeasure`\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.modifyprobability!",
    "page": "Reference",
    "title": "SDDP.modifyprobability!",
    "category": "Function",
    "text": "modifyprobability!(measure::AbstractRiskMeasure,\n        riskadjusted_distribution,\n        original_distribution::Vector{Float64},\n        observations::Vector{Float64},\n        m::SDDPModel,\n        sp::JuMP.Model\n)\n\nDescription\n\nCalculate the risk-adjusted probability of each scenario using the 'change-of-probabilities' approach of Philpott, de Matos, and Finardi,(2013). On solving multistage stochastic programs with coherent risk measures. Operations Research 61(4), 957-970.\n\nArguments\n\nmeasure::AbstractRiskMeasure\n\nThe risk measure\n\nriskadjusted_distribution\n\nA new probability distribution\n\noriginal_distribution::Vector{Float64}\n\nThe original probability distribution.\n\nobservations::Vector{Float64}\n\nThe vector of objective values from the next stage  problems (one for each scenario).\n\nm::SDDPModel\n\nThe full SDDP model\n\nsp::JuMP.Model\n\nThe stage problem that the cut will be added to.\n\n\n\n"
},

{
    "location": "apireference.html#Risk-Measures-1",
    "page": "Reference",
    "title": "Risk Measures",
    "category": "section",
    "text": "AbstractRiskMeasure\nExpectation\nNestedAVaR\nSDDP.modifyprobability!"
},

{
    "location": "apireference.html#Cut-Oracles-1",
    "page": "Reference",
    "title": "Cut Oracles",
    "category": "section",
    "text": "AbstractCutOracle\nDefaultCutOracle\nLevelOneCutOracle\nSDDP.storecut!!\nSDDP.validcuts"
},

{
    "location": "apireference.html#Solve-1",
    "page": "Reference",
    "title": "Solve",
    "category": "section",
    "text": "solve\nMonteCarloSimulation\nBoundConvergence\nSerial\nAsyncronous"
},

{
    "location": "apireference.html#Results-1",
    "page": "Reference",
    "title": "Results",
    "category": "section",
    "text": "getbound\nsimulate\n@visualise"
},

{
    "location": "apireference.html#Save/Load-1",
    "page": "Reference",
    "title": "Save/Load",
    "category": "section",
    "text": "loadcuts!"
},

]}
