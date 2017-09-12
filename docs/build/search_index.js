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
    "location": "quick.html#Initialising-the-model-object-1",
    "page": "Quick Start",
    "title": "Initialising the model object",
    "category": "section",
    "text": "The first step is to initialise the SDDP model object. We do this using the following syntax:If we have more than one markov state:m = SDDPModel([;kwargs...]) do sp, stage, markov_state\n\n    # Stage problem definition where `sp` is a `JuMP.Model`object,\n\nendOtherwise if we have a single markov statem = SDDPModel([;kwargs...]) do sp, stage\n\n  # Stage problem definition\n\nend"
},

{
    "location": "quick.html#Solve-1",
    "page": "Quick Start",
    "title": "Solve",
    "category": "section",
    "text": "To solve the SDDP model m we use status = solve(m::SDDPModel; kwargs...). This accepts a number of keyword arguments to control the solution process."
},

{
    "location": "quick.html#Simulate-1",
    "page": "Quick Start",
    "title": "Simulate",
    "category": "section",
    "text": ""
},

{
    "location": "quick.html#Visualise-1",
    "page": "Quick Start",
    "title": "Visualise",
    "category": "section",
    "text": ""
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
    "location": "apireference.html#SDDP.SDDPModel",
    "page": "Reference",
    "title": "SDDP.SDDPModel",
    "category": "Type",
    "text": "SDDPModel(;kwargs...) do ...\n\nend\n\nDescription\n\nThis function constructs an SDDPModel. SDDPModel takes the following keyword arguments. Some are required, and some are optional.\n\nRequired Keyword arguments\n\nstages::Int\n\nThe number of stages in the problem. A stage is defined as each step in time at  which a decion can be made. Defaults to 1.\n\nobjective_bound::Float64\n\nA valid bound on the initial value/cost to go. i.e. for maximisation this may  be some large positive number, for minimisation this may be some large negative  number.\n\nsolver::MathProgBase.AbstractMathProgSolver\n\nMathProgBase compliant solver that returns duals from a linear program. If this  isn't specified then you must use JuMP.setsolver(sp, solver) in the stage  definition.\n\nOptional Keyword arguments\n\nsense\n\nMust be either :Max or :Min. Defaults to :Min.\n\ncut_oracle::SDDP.AbstractCutOracle\n\nThe cut oracle is responsible for collecting and storing the cuts that define  a value function. The cut oracle may decide that only a subset of the total  discovered cuts are relevant, which improves solution speed by reducing the  size of the subproblems that need solving. Currently must be one of     * DefaultCutOracle() (see DefaultCutOracle for explanation)     * LevelOneCutOracle()(see LevelOneCutOracle for explanation)\n\nrisk_measure\n\nIf a single risk measure is given (i.e. risk_measure = Expectation()), then  this measure will be applied to every stage in the problem. Another option is  to provide a vector of risk measures. There must be one element for every  stage. For example:\n\nrisk_measure = [ NestedAVaR(lambda=0.5, beta=0.25), Expectation() ]\n\nwill apply the i'th element of risk_measure to every Markov state in the i'th stage. The last option is to provide a vector (one element for each stage) of vectors of risk measures (one for each Markov state in the stage). For example:\n\nrisk_measure = [\n# Stage 1 Markov 1 # Stage 1 Markov 2 #\n    [ Expectation(), Expectation() ],\n    # ------- Stage 2 Markov 1 ------- ## ------- Stage 2 Markov 2 ------- #\n    [ NestedAVaR(lambda=0.5, beta=0.25), NestedAVaR(lambda=0.25, beta=0.3) ]\n    ]\n\nNote that even though the last stage does not have a future cost function associated with it (as it has no children), we still have to specify a risk measure. This is necessary to simplify the implementation of the algorithm.\n\nFor more help see NestedAVaR or Expectation.\n\nmarkov_transition\n\nDefine the transition probabilties of the stage graph. If a single array is given, it is assumed that there is an equal number of Markov states in each stage and the transition probabilities are stage invariant. Row indices represent the Markov state in the previous stage. Column indices represent the Markov state in the current stage. Therefore:\n\nmarkov_transition = [0.1 0.9; 0.8 0.2]\n\nis the transition matrix when there is 10% chance of transitioning from Markov state 1 to Markov state 1, a 90% chance of transitioning from Markov state 1 to Markov state 2, an 80% chance of transitioning from Markov state 2 to Markov state 1, and a 20% chance of transitioning from Markov state 2 to Markov state 2.\n\nReturns\n\nm: the SDDPModel\n\n\n\n"
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
    "location": "apireference.html#SDDP.@stageobjective",
    "page": "Reference",
    "title": "SDDP.@stageobjective",
    "category": "Macro",
    "text": "@stageobjective!(sp, kw=noises, objective)\n\nDescription\n\nDefine an objective that depends on the realization of the stagewise noise. objective can be any valid third argument to the JuMP @objective macro (i.e. @objective(sp, Min, objective)) that utilises the variable kw that takes the realizations defined in noises.\n\nExamples\n\n@stageobjective(sp, w=1:2, w * x)\n@stageobjective(sp, i=1:2, w[i]^2 * x)\n@stageobjective(sp, i=1:2, x[i])\n\n\n\n"
},

{
    "location": "apireference.html#Objective-1",
    "page": "Reference",
    "title": "Objective",
    "category": "section",
    "text": "@stageobjective"
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
    "text": "NestedAVaR(;lambda=1.0, beta=1.0)\n\nDescription\n\nA risk measure that is a convex combination of Expectation and Average Value @ Risk (also called Conditional Value @ Risk).\n\nλ * E[x] + (1 - λ) * AV@R(1-β)[x]\n\nKeyword Arguments\n\nlambda\n\nConvex weight on the expectation ((1-lambda) weight is put on the AV@R component. Inreasing values of lambda are less risk averse (more weight on expecattion)\n\nbeta\n\nThe quantile at which to calculate the Average Value @ Risk. Increasing values  of beta are less risk averse. If beta=0, then the AV@R component is the  worst case risk measure.\n\nReturns\n\nm::NestedAVaR<:AbstractRiskMeasure\n\n\n\n"
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
    "text": "AbstractRiskMeasure\nExpectation\nNestedAVaR\nmodifyprobability!"
},

{
    "location": "apireference.html#SDDP.AbstractCutOracle",
    "page": "Reference",
    "title": "SDDP.AbstractCutOracle",
    "category": "Type",
    "text": "AbstractCutOracle\n\nDescription\n\nAbstract type for all cut oracles.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.DefaultCutOracle",
    "page": "Reference",
    "title": "SDDP.DefaultCutOracle",
    "category": "Type",
    "text": "DefaultCutOracle()\n\nDescription\n\nInitialize the default cut oracle.\n\nThis oracle keeps every cut discovered and does not perform cut selection.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.LevelOneCutOracle",
    "page": "Reference",
    "title": "SDDP.LevelOneCutOracle",
    "category": "Function",
    "text": "LevelOneCutOracle()\n\nDescription\n\nInitialize the cut oracle for Level One cut selection. See:\n\nV. de Matos,A. Philpott, E. Finardi, Improving the performance of Stochastic Dual Dynamic Programming, Journal of Computational and Applied Mathematics 290 (2015) 196–208.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.storecut!",
    "page": "Reference",
    "title": "SDDP.storecut!",
    "category": "Function",
    "text": "storecut!(oracle::AbstactCutOracle, m::SDDPModel, sp::JuMP.Model, cut::Cut)\n\nDescription\n\nStore the cut cut in the Cut Oracle oracle. oracle will belong to the subproblem sp in the SDDPModel m.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.validcuts",
    "page": "Reference",
    "title": "SDDP.validcuts",
    "category": "Function",
    "text": "validcuts(oracle::AbstactCutOracle)\n\nDescription\n\nReturn an iterable list of all the valid cuts contained within oracle.\n\n\n\n"
},

{
    "location": "apireference.html#Cut-Oracles-1",
    "page": "Reference",
    "title": "Cut Oracles",
    "category": "section",
    "text": "AbstractCutOracle\nDefaultCutOracle\nLevelOneCutOracle\nstorecut!\nvalidcuts"
},

{
    "location": "apireference.html#JuMP.solve",
    "page": "Reference",
    "title": "JuMP.solve",
    "category": "Function",
    "text": "solve(m::SDDPModel; kwargs...)\n\nDescription\n\nSolve the SDDPModel m using SDDP. Accepts a number of keyword arguments to control the solution process.\n\nPositional arguments\n\nm: the SDDPModel to solve\n\nKeyword arguments\n\nmax_iterations::Int:  The maximum number of cuts to add to a single stage problem before terminating.  Defaults to 10.\ntime_limit::Real:  The maximum number of seconds (in real time) to compute for before termination.  Defaults to Inf.\nsimulation::MonteCarloSimulation: see MonteCarloSimulation\nbound_convergence::BoundConvergence: see BoundConvergence\ncut_selection_frequency::Int:  Frequency (by iteration) with which to rebuild subproblems using a subset of  cuts. Frequent cut selection (i.e. cut_selection_frequency is small) reduces  the size of the subproblems that are solved, but incurrs the overhead of rebuilding  the subproblems. However, infrequent cut selection (i.e.  cut_selection_frequency is large) allows the subproblems to grow large (many  constraints) leading to an increase in the solve time of individual subproblems.  Defaults to 0 (never run).\nprint_level::Int:   0 - off: nothing logged to screen (still written to log file if specified).   1 - on: solve iterations written to screen.   Defaults to 1\nlog_file::String:  Relative filename to write the log to disk. Defaults to \"\" (no log written)\nsolve_type:  One of\nAsyncronous() - solve using a parallelised algorithm\nSerial() - solve using a serial algorithm\nDefault chooses automatically based on the number of available processors.\nreduce_memory_footprint::Bool:  Implements the idea proposed in https://github.com/JuliaOpt/JuMP.jl/issues/969#issuecomment-282191105  to reduce the memory consumption when running SDDP. This is an issue if you  wish to save the model m to disk since it discards important information.  Defaults to false.\ncut_output_file::String:  Relative filename to write discovered cuts to disk. Defaults to \"\" (no cuts written)\n\nReturns\n\nstatus::Symbol:  Reason for termination. One of\n:solving\n:interrupted\n:converged\n:max_iterations\n:bound_convergence\n:time_limit\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.MonteCarloSimulation",
    "page": "Reference",
    "title": "SDDP.MonteCarloSimulation",
    "category": "Type",
    "text": "MonteCarloSimulation(;kwargs...)\n\nDescription\n\nCollection of settings to control the simulation phase of the SDDP solution process.\n\nArguments\n\nfrequency::Int\n\nThe frequency (by iteration) with which to run the policy simulation phase of the algorithm in order to construct a statistical bound for the policy. Defaults to 0 (never run).\n\nmin::Float64\n\nMinimum number of simulations to conduct before constructing a confidence interval for the bound. Defaults to 20.\n\nstep::Float64\n\nNumber of additional simulations to conduct before constructing a new confidence interval for the bound. Defaults to 1.\n\nmax::Float64\n\nMaximum number of simulations to conduct in the policy simulation phase. Defaults to min.\n\nconfidence::Float64\n\nConfidence level of the confidence interval. Defaults to 0.95 (95% CI).\n\ntermination::Bool\n\nWhether to terminate the solution algorithm with the status :converged if the deterministic bound is with in the statistical bound after max simulations. Defaults to false.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.BoundConvergence",
    "page": "Reference",
    "title": "SDDP.BoundConvergence",
    "category": "Type",
    "text": "BoundConvergence(;kwargs...)\n\nDescription\n\nCollection of settings to control the bound stalling convergence test.\n\nArguments\n\niterations::Int\n\nTerminate if the maximum deviation in the deterministic bound from the mean over the last iterations number of iterations is less than rtol (in relative terms) or atol (in absolute terms).\n\nrtol::Float64\n\nMaximum allowed relative deviation from the mean. Defaults to 0.0\n\natol::Float64\n\nMaximum allowed absolute deviation from the mean. Defaults to 0.0\n\n\n\n"
},

{
    "location": "apireference.html#Solve-1",
    "page": "Reference",
    "title": "Solve",
    "category": "section",
    "text": "solve\nMonteCarloSimulation\nBoundConvergence"
},

{
    "location": "apireference.html#SDDP.simulate",
    "page": "Reference",
    "title": "SDDP.simulate",
    "category": "Function",
    "text": "simulate(m::SDDPPModel,variables::Vector{Symbol};\n    noises::Vector{Int}, markovstates::Vector{Int})\n\nDescription\n\nPerform a historical simulation of the current policy in model  m.\n\nnoises is a vector with one element for each stage giving the index of the (in-sample) stagewise independent noise to sample in each stage. markovstates is a vector with one element for each stage giving the index of the (in-sample) markov state to sample in each stage.\n\nExamples\n\nsimulate(m, [:x, :u], noises=[1,2,2], markovstates=[1,1,2])\n\n\n\nresults = simulate(m::SDDPPModel, N::Int, variables::Vector{Symbol})\n\nDescription\n\nPerform N Monte-Carlo simulations of the current policy in model m saving the values of the variables named in variables at every stage.\n\nresults is a vector containing a dictionary for each simulation. In addition to the variables specified in the function call, other special keys are:\n\n:stageobjective - costs incurred during the stage (not future)\n:obj            - objective of the stage including future cost\n:markov         - index of markov state visited\n:noise          - index of noise visited\n:objective      - Total objective of simulation\n\nAll values can be accessed as follows\n\nresults[simulation index][key][stage]\n\nwith the exception of :objective which is just\n\nresults[simulation index][:objective]\n\nExamples\n\nresults = simulate(m, 10, [:x, :u])\nresults[1][:objective] # objective of simulation 1\nmean(r[:objective] for r in results) # mean objective of the simulations\nresults[2][:x][3] # value of :x in stage 3 in second simulation\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.getbound",
    "page": "Reference",
    "title": "SDDP.getbound",
    "category": "Function",
    "text": "getbound(m)\n\nDescription\n\nGet the lower (if minimizing), or upper (if maximizing) bound of the solved SDDP model m.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.newplot",
    "page": "Reference",
    "title": "SDDP.newplot",
    "category": "Function",
    "text": "SDDP.newplot()\n\nDescription\n\nInitialize a new SimulationPlot.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.addplot!",
    "page": "Reference",
    "title": "SDDP.addplot!",
    "category": "Function",
    "text": "SDDP.addplot!(p::SimulationPlot, ivals::AbstractVector{Int}, tvals::AbstractVector{Int}, f::Function; kwargs...)\n\nDescription\n\nAdd a new figure to the SimulationPlot p, where the y-value is given by f(i, t) for all i in ivals (one for each series) and t in tvals (one for each stage).\n\nKeywords\n\nxlabel: set the xaxis label;\nylabel: set the yaxis label;\ntitle: set the title of the plot;\nymin: set the minimum y value;\nymax: set the maximum y value;\ncumulative: plot the additive accumulation of the value across the stages;\ninterpolate: interpolation method for lines between stages. Defaults to \"linear\"  see the d3 docs\n\nfor all options.\n\nExamples\n\nresults = simulate(m, 10)\np = SDDP.newplot()\nSDDP.addplot!(p, 1:10, 1:3, (i,t)->results[i][:stageobjective][t])\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.show",
    "page": "Reference",
    "title": "SDDP.show",
    "category": "Function",
    "text": "show(p::SimulationPlot)\n\nDescription\n\nLaunch a browser and render the SimulationPlot plot p.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.plotvaluefunction",
    "page": "Reference",
    "title": "SDDP.plotvaluefunction",
    "category": "Function",
    "text": " SDDP.plotvaluefunction(m::SDDPModel, stage::Int, markovstate::Int, states::Union{Float64, AbstractVector{Float64}}...; label1=\"State 1\", label2=\"State 2\")\n\nDescription\n\nPlot the value function of stage stage and Markov state markovstate in the SDDPModel m at the points in the discretized state space given by states. If the value in states is a real number, the state is evaluated at that point. If the value is a vector, the state is evaluated at all the points in the vector. At most two states can be vectors.\n\nExamples\n\nSDDP.plotvaluefunction(m, 2, 1, 0.0:0.1:1.0, 0.5, 0.0:0.1:1.0; label1=\"State 1\", label2=\"State 3\")\n\n\n\n SDDP.plotvaluefunction(vf::DefaultValueFunction, is_minimization::Bool, states::Union{Float64, AbstractVector{Float64}}...; label1=\"State 1\", label2=\"State 2\")\n\nDescription\n\nPlot the value function vf at the points in the discretized state space given by states. If the value in states is a real number, the state is evaluated at that point. If the value is a vector, the state is evaluated at all the points in the vector. At most two states can be vectors. is_minimization is true if the subproblem is a minimization and false otherwise.\n\nExamples\n\nSDDP.plotvaluefunction(vf, 0.0:0.1:1.0, 0.5, 0.0:0.1:1.0; label1=\"State 1\", label2=\"State 3\")\n\n\n\n"
},

{
    "location": "apireference.html#Results-1",
    "page": "Reference",
    "title": "Results",
    "category": "section",
    "text": "simulate\ngetbound\nnewplot\naddplot!\nshow\nplotvaluefunction"
},

{
    "location": "apireference.html#Save/Load-1",
    "page": "Reference",
    "title": "Save/Load",
    "category": "section",
    "text": "<!– loadcuts! –>"
},

]}
