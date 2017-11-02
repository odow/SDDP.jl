var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Manual",
    "title": "Manual",
    "category": "page",
    "text": "CurrentModule = SDDP"
},

{
    "location": "index.html#Overview-1",
    "page": "Manual",
    "title": "Overview",
    "category": "section",
    "text": "SDDP.jl - Stochastic Dual Dynamic Programming in Julia.SDDP.jl is a package for solving large multi-stage convex stocastic optimization problems. In this manual, we're going to assume a reasonable amount of background knowledge about stochastic optimization, the SDDP algorithm, Julia, and JuMP.note: Note\nIf you don't have that background, you may want to brush up on some reading material."
},

{
    "location": "index.html#Types-of-problems-SDDP.jl-can-solve-1",
    "page": "Manual",
    "title": "Types of problems SDDP.jl can solve",
    "category": "section",
    "text": "To start, lets discuss the types of problems SDDP.jl can solve, since it has a few features that are non-standard, and it is missing some features that are standard.SDDP.jl can solve multi-stage convex stochastic optimizations problems witha finite discrete number of states;\ncontinuous state and control variables;\nHazard-Decision (Wait-and-See) uncertainty realization;\nstagewise independent uncertainty in the RHS of the constraints that is  drawn from a finite discrete distribution;\nstagewise independent uncertainty in the objective function that is  drawn from a finite discrete distribution;\na markov chain for temporal dependence. The markov chain forms a directed,  acyclic, feed-forward graph with a finite (and at least one) number of  markov states in each stage.note: Note\nStagewise independent uncertainty in the constraint coefficients is not supported. You should reformulate the problem, or model the uncertainty as a markov chain.In this manual, we detail the many features of SDDP.jl through the classic example of balancing a portfolio of stocks and bonds over time."
},

{
    "location": "index.html#Getting-started-1",
    "page": "Manual",
    "title": "Getting started",
    "category": "section",
    "text": "This package is unregistered so you will need to Pkg.clone it as follows:Pkg.clone(\"https://github.com/odow/SDDP.jl.git\")If you want to use the parallel features of SDDP.jl, you should start Julia with some worker processes (julia -p N), or add by running julia> addprocs(N) in a running Julia session."
},

{
    "location": "index.html#Formulating-the-problem-1",
    "page": "Manual",
    "title": "Formulating the problem",
    "category": "section",
    "text": "... still to do ...For now, go look at the examples."
},

{
    "location": "index.html#The-Asset-Management-Problem-1",
    "page": "Manual",
    "title": "The Asset Management Problem",
    "category": "section",
    "text": "The goal of the asset management problem is to choose an investment portfolio that is composed of stocks and bonds in order to meet a target wealth goal at the end of the time horizon. After five, and ten years, the agent observes the portfolio and is able to re-balance their wealth between the two asset classes. As an extension to the original problem, we introduce two new random variables. The first that represents a source of additional wealth in years 5 and 10. The second is an immediate reward that the agent incurs for holding stocks at the end of years 5 and 10. This can be though of as a dividend that cannot be reinvested."
},

{
    "location": "index.html#Communicating-the-problem-to-the-solver-1",
    "page": "Manual",
    "title": "Communicating the problem to the solver",
    "category": "section",
    "text": "The second step in the optimization process is communicating the problem to the solver. To do this, we are going to build each subproblem as a JuMP model, and provide some metadata that describes how the JuMP subproblems inter-relate."
},

{
    "location": "index.html#The-Model-Constructor-1",
    "page": "Manual",
    "title": "The Model Constructor",
    "category": "section",
    "text": "The core type of SDDP.jl is the SDDPModel object. It can be constructed withm = SDDPModel( ... metadata as keyword arguments ... ) do sp, t, i\n    ... JuMP subproblem definition ...\nendWe draw the readers attention to three sections in the SDDPModel constructor."
},

{
    "location": "index.html#Keyword-Metadata-1",
    "page": "Manual",
    "title": "Keyword Metadata",
    "category": "section",
    "text": "For a comprehensive list of options, checkout the SDDPModel API or type julia> ? SDDPModel into a Julia REPL. However, we'll briefly list the important ones here.Required Keyword argumentsstages::Int: the number of stages in the problem. A stage is defined as  each step in time at which a decion can be made. Defaults to 1.objective_bound: a valid bound on the initial value/cost to go. i.e. for maximisation this may be some large positive number, for minimisation this may be some large negative number. You can ether pass a single value (used for all stages), a vector of values (one for each stage), or a vector of vectors of values (one vector for each stage, containing a vector with one element for each markov stage). Another option is to pass a function that takes two inputs so that f(t, i) returns the valid bound for stage t and markov state i.solver::MathProgBase.AbstractMathProgSolver: any MathProgBase compliant solver that returns duals from a linear program. If this isn't specified then you must use JuMP.setsolver(sp, solver) in the stage definition.Optional Keyword argumentssense: must be either :Max or :Min. Defaults to :Min.cut_oracle: the cut oracle is responsible for collecting and storing the cuts that define a value function. The cut oracle may decide that only a subset of the total discovered cuts are relevant, which improves solution speed by reducing the size of the subproblems that need solving. Currently must be one of\nDefaultCutOracle() (see DefaultCutOracle for explanation)\nLevelOneCutOracle()(see LevelOneCutOracle for explanation)risk_measure: if a single risk measure is given (i.e.  risk_measure = Expectation()), then this measure will be applied to every  stage in the problem. Another option is to provide a vector of risk  measures. There must be one element for every stage. For example:risk_measure = [ NestedAVaR(lambda=0.5, beta=0.25), Expectation() ]will apply the i'th element of risk_measure to every Markov state in the    i'th stage. The last option is to provide a vector (one element for each    stage) of vectors of risk measures (one for each Markov state in the stage).    For example:risk_measure = [\n# Stage 1 Markov 1 # Stage 1 Markov 2 #\n   [ Expectation(), Expectation() ],\n   # ------- Stage 2 Markov 1 ------- ## ------- Stage 2 Markov 2 ------- #\n   [ NestedAVaR(lambda=0.5, beta=0.25), NestedAVaR(lambda=0.25, beta=0.3) ]\n   ]Like the objective bound, another option is to pass a function that takes\nas arguments the stage and markov state and returns a risk measure.Note that even though the last stage does not have a future cost function    associated with it (as it has no children), we still have to specify a risk    measure. This is necessary to simplify the implementation of the algorithm.For more help see NestedAVaR or Expectation.markov_transition: define the transition probabilties of the stage graph. If a single array is given, it is assumed that there is an equal number of Markov states in each stage and the transition probabilities are stage invariant. Row indices represent the Markov state in the previous stage. Column indices represent the Markov state in the current stage. Therefore:markov_transition = [0.1 0.9; 0.8 0.2]is the transition matrix when there is 10% chance of transitioning from Markov    state 1 to Markov state 1, a 90% chance of transitioning from Markov state 1    to Markov state 2, an 80% chance of transitioning from Markov state 2 to Markov    state 1, and a 20% chance of transitioning from Markov state 2 to Markov state 2."
},

{
    "location": "index.html#do-sp,-t,-i-...-end-1",
    "page": "Manual",
    "title": "do sp, t, i ... end",
    "category": "section",
    "text": "This constructor is just syntactic sugar to make the process of defining a model a little tidier. It's nothing special to SDDP.jl and many users will be familiar with it (for example, the open(file, \"w\") do io ... end syntax for file IO).An anonymous function with three arguments (sp, t and i, although these can be named arbitrarily) is constructed. The body of the function should build the subproblem as the JuMP model sp for stage t and markov state i. t is an integer that ranges from 1 to the number of stages. i is also an integer that ranges from 1 to the number of markov states in stage t.Users are also free to explicitly construct a function that takes three arguments, and pass that as the first argument to SDDPModel, along with the keyword arguments. For example:function buildsubproblem!(sp::JuMP.Model, t::Int, i::Int)\n    ... define states etc. ...\nend\nm = SDDPModel(buildsubproblem!; ... metadata as keyword arguments ...)note: Note\nIf you don't have any markov states in the model, you don't have to include the third argument in the constructor. SDDPModel() do sp, t ... end is also valid syntax."
},

{
    "location": "index.html#JuMP-Subproblem-1",
    "page": "Manual",
    "title": "JuMP Subproblem",
    "category": "section",
    "text": "In the next sections, we explain in detail how for model state variables, constraints, the stage objective, and any uncertainties in the model. However, you should keep a few things in mind:the body of the do sp, t, i ... end block is just a normal Julia function body. As such, standard scoping rules apply.\nyou can use t and i whenever, and however, you like. For example:m = SDDPModel() do sp, t, i\n    if t == 1\n        # do something in the first stage only\n    else\n        # do something in the other stages\n    end\nendsp is just a normal JuMP model. You could (if so desired), set the solve hook, or add quadratic constraints (provided you have a quadratic solver)."
},

{
    "location": "index.html#State-Variables-1",
    "page": "Manual",
    "title": "State Variables",
    "category": "section",
    "text": "We can define a new state variable in the stage problem sp using the @state macro:@state(sp, x >= 0.5, x0==1)The second argument (x) refers to the outgoing state variable (i.e. the value at the end of the stage). The third argument (x0) refers to the incoming state variable (i.e. the value at the beginning of the stage). For users familiar with SDDP, SDDP.jl handles all the calculation of the dual variables needed to evaluate the cuts automatically behind the scenes.The @state macro is just short-hand for writing:@variable(sp, x >= 0.5)\n@variable(sp, x0, start=1)\nSDDP.statevariable!(sp, x0, x)note: Note\nThe start=1 is only every used behind the scenes by the first stage problem. It's really just a nice syntatic trick we use to make specifying the model a bit more compact.This illustrates how we can use indexing just as we would in a JuMP @variable macro:X0 = [3.0, 2.0]\n@state(sp, x[i=1:2], x0==X0[i])In this case, both x and x0 are JuMP dicts that can be indexed with the keys 1 and 2. All the indices must be specified in the second argument, but they can be referred to in the third argument. The indexing of x0 will be identical to that of x.There is also a plural version of the @state macro:@states(sp, begin\n    x >= 0.0, x0==1\n    y >= 0.0, y0==1\nend)"
},

{
    "location": "index.html#Standard-JuMP-machinery-1",
    "page": "Manual",
    "title": "Standard JuMP machinery",
    "category": "section",
    "text": "Remember that sp is just a normal JuMP model, and so (almost) anything you can do in JuMP, you can do in SDDP.jl. The one exception is the objective, which we detail in the next section.However, control variables are just normal JuMP variables and can be created using @variable or @variables. Dynamical constraints, and feasiblity sets can be specified using @constraint or @constraints."
},

{
    "location": "index.html#The-stage-objective-1",
    "page": "Manual",
    "title": "The stage objective",
    "category": "section",
    "text": "If there is no stagewise independent uncertainty in the objective, then the stage objective (i.e. ignoring the future cost) can be set via the @stageobjective macro. This is similar to the JuMP @objective macro, but without the sense argument. For example:@stageobjective(sp, obj)If there is stagewise independent noise in the objective, we add an additional argument to @stageobjective that has the form kw=realizations.kw is a symbol that can appear anywhere in obj, and realizations is a vector of realizations of the uncertainty. For example:@stageobjective(sp, kw=realizations, obj)\nsetnoiseprobability!(sp, [0.2, 0.3, 0.5])setnoiseprobability! can be used to specify the finite discrete distribution of the realizations (it must sum to 1.0). If you don't explicitly call setnoiseprobability!, the distribution is assumed to be uniform.Other examples include:# noise is a coefficient\n@stageobjective(sp, c=[1.0, 2.0, 3.0], c * x)\n# noise is used to index a variable\n@stageobjective(sp, i=[1,2,3], 2 * x[i])"
},

{
    "location": "index.html#Dynamics-with-linear-noise-1",
    "page": "Manual",
    "title": "Dynamics with linear noise",
    "category": "section",
    "text": "SDDP.jl also supports uncertainty in the right-hand-side of constraints. Instead of using the JuMP @constraint macro, we need to use the @rhsnoise macro:@rhsnoise(sp, w=[1,2,3], x <= w)\nsetnoiseprobability!(sp, [0.2, 0.3, 0.5])Compared to @constraint, there are a couple of notable differences:indexing is not supported;\nthe second argument is a kw=realizations key-value pair like the @stageobjective;\nthe kw can on either side of the constraint as written, but when normalised  to an Ax <= b form, it must only appear in the b vector.Multiple @rhsnoise constraints can be added, however they must have an identical number of elements in the realizations vector.For example, the following are invalid in SDDP:# noise appears as a variable coefficient\n@rhsnoise(sp, w=[1,2,3], w * x <= 1)\n\n# JuMP style indexing\n@rhsnoise(sp, w=[1,2,3], [i=1:10; mod(i, 2) == 0], x[i] <= w)\n\n# noises have different number of realizations\n@rhsnoise(sp, w=[1,2,3], x <= w)\n@rhsnoise(sp, w=[2,3],   x >= w-1)note: Note\nNoises in the constraints are sampled with the noise in the objective. Therefore, there should be the same number of elements in the realizations for the stage objective, as there are in the constraint noise.There is also a plural form of the @rhsnoise macro:@rhsnoises(sp, w=[1,2,3], begin\n    x <= w\n    x >= w-1\nend)\nsetnoiseprobability!(sp, [0.2, 0.3, 0.5])"
},

{
    "location": "index.html#Asset-Management-Example-1",
    "page": "Manual",
    "title": "Asset Management Example",
    "category": "section",
    "text": "We now have all the information necessary to define the Asset Management example in SDDP.jl:using SDDP, JuMP, Clp\n\nm = SDDPModel(\n               # we are minimizing\n                sense = :Min,\n               # there are 4 stages\n               stages = 4,\n               # a large negative value\n      objective_bound = -1000.0,\n               # a MathOptBase compatible solver\n               solver = ClpSolver(),\n               # transition probabilities of the lattice\n    markov_transition = Array{Float64, 2}[\n                        [1.0]',\n                        [0.5 0.5],\n                        [0.5 0.5; 0.5 0.5],\n                        [0.5 0.5; 0.5 0.5]\n                    ],\n               # risk measures for each stage\n         risk_measure = [\n                        Expectation(),\n                        Expectation(),\n                        NestedAVaR(lambda = 0.5, beta=0.5),\n                        Expectation()\n                    ]\n                            ) do sp, t, i\n    # Data used in the problem\n    ωs = [1.25, 1.06]\n    ωb = [1.14, 1.12]\n    Φ  = [-1, 5]\n    Ψ  = [0.02, 0.0]\n\n    # state variables\n    @states(sp, begin\n        xs >= 0, xsbar==0\n        xb >= 0, xbbar==0\n    end)\n\n    if t == 1 # we can switch based on the stage\n        # a constraint without noise\n        @constraint(sp, xs + xb == 55 + xsbar + xbbar)\n        # an objective without noise\n        @stageobjective(sp, 0)\n    elseif t == 2 || t == 3\n        # a constraint with noisein the RHS\n        @rhsnoise(sp, φ=Φ, ωs[i] * xsbar + ωb[i] * xbbar + φ == xs + xb)\n        # an objective with noise\n        @stageobjective(sp, ψ = Ψ, -ψ * xs)\n        # noise will be sampled as (Φ[1], Ψ[1]) w.p. 0.6, (Φ[2], Ψ[2]) w.p. 0.4\n        setnoiseprobability!(sp, [0.6, 0.4])\n    else # when t == 4\n        # some control variables\n        @variable(sp, u >= 0)\n        @variable(sp, v >= 0)\n        # dynamics constraint\n        @constraint(sp, ωs[i] * xsbar + ωb[i] * xbbar + u - v == 80)\n        # an objective without noise\n        @stageobjective(sp, 4u - v)\n    end\nend"
},

{
    "location": "index.html#Solving-the-problem-efficiently-1",
    "page": "Manual",
    "title": "Solving the problem efficiently",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Understanding-the-solution-1",
    "page": "Manual",
    "title": "Understanding the solution",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Simulating-the-policy-1",
    "page": "Manual",
    "title": "Simulating the policy",
    "category": "section",
    "text": "You can perform a Monte-Carlo simulation of the policy using the simulate function:simulationresults = simulate(m, 100, [:xs, :xb])simulationresults is a vector of dictionaries (one for each simulation). It has, as keys, a vector (with one element for each stage) of the optimal solution of xs and xb. For example, to query the value of xs in the third stage of the tenth simulation, we can call:simulationresults[10][:xs][3]Alternatively, you can peform one simulation with a given realization for the Markov and stagewise independent noises:simulationresults = simulate(m, [:xs], markov=[1,2,1,1], noise=[1,2,2,3])"
},

{
    "location": "index.html#Visualizing-the-policy-simulation-1",
    "page": "Manual",
    "title": "Visualizing the policy simulation",
    "category": "section",
    "text": "First, we create a new plotting object with SDDP.newplot(). Next, we can add any number of subplots to the visualization via the SDDP.addplot! function. Finally, we can launch a web browser to display the plot with SDDP.show.See the SDDP.addplot! documentation for more detail."
},

{
    "location": "index.html#Visualizing-the-Value-Function-1",
    "page": "Manual",
    "title": "Visualizing the Value Function",
    "category": "section",
    "text": "Another way to understand the solution is to project the value function into 3 dimensions. This can be done using the method SDDP.plotvaluefunction.SDDP.plotvaluefunction(m, 1, 1, 0:1.0:100, 0:1.0:100;\n    label1=\"Stocks\", label2=\"Bonds\")This will open up a web browser and display a Plotly figure that looks similar to (Image: 3-Dimensional visualisation of Value Function)"
},

{
    "location": "index.html#Saving-models-1",
    "page": "Manual",
    "title": "Saving models",
    "category": "section",
    "text": "Saving a model is as simple as calling:SDDP.savemodel!(\"<filename>\", m)Later, you can run:m = SDDP.loadmodel(\"<filename>\")note: Note\nSDDP.savemodel! relies on the base Julia serialize function. This is not backwards compatible with previous versions of Julia, or guaranteed to be forward compatible with future versions. You should only use this to save models for short periods of time. Don't save a model you want to come back to in a year.Another (more persistent) method is to use the cut_output_file keyword option in SDDP.solve. This will create a csv file containing a list of all the cuts. These can be loaded at a later date usingSDDP.loadcuts!(m, \"<filename>\")"
},

{
    "location": "index.html#Extras-for-experts-1",
    "page": "Manual",
    "title": "Extras for experts",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#New-risk-measures-1",
    "page": "Manual",
    "title": "New risk measures",
    "category": "section",
    "text": "SDDP.jl makes it easy to create new risk measures. First, create a new subtype of the abstract type SDDP.AbstractRiskMeasure:immutable MyNewRiskMeasure <: SDDP.AbstractRiskMeasure\nendThen, overload the method SDDP.modifyprobability! for your new type. SDDP.modifyprobability! has the following signature:SDDP.modifyprobability!(\n        measure::AbstractRiskMeasure,\n        riskadjusted_distribution,\n        original_distribution::Vector{Float64},\n        observations::Vector{Float64},\n        m::SDDPModel,\n        sp::JuMP.Model\n)where original_distribution contains the risk netural probability of each outcome in observations occurring (so that the probability of observations[i] is original_distribution[i]). The method should modify (in-place) the elements of riskadjusted_distribution to represent the risk-adjusted probabilities of the distribution.To illustrate this, we shall define the worst-case riskmeasure (which places all the probability on the worst outcome):immutable WorstCase <: SDDP.AbstractRiskMeasure end\nfunction SDDP.modifyprobability!(measure::WorstCase,\n        riskadjusted_distribution,\n        original_distribution::Vector{Float64},\n        observations::Vector{Float64},\n        m::SDDPModel,\n        sp::JuMP.Model\n    )\n    if JuMP.getobjectivesense(sp) == :Min\n        # if minimizing, worst is the maximum outcome\n        idx = indmax(observations)\n    else\n        # if maximizing, worst is the minimum outcome\n        idx = indmin(observations)\n    end\n    # zero all probabilities\n    riskadjusted_distribution .= 0.0\n    # set worst to 1.0\n    riskadjusted_distribution[idx] = 1.0\n    # return\n    return nothing\nendThe risk measure WorstCase() can now be used in any SDDP model.note: Note\nThis method gets called a lot, so the usual Julia performance tips apply."
},

{
    "location": "index.html#New-cut-oracles-1",
    "page": "Manual",
    "title": "New cut oracles",
    "category": "section",
    "text": ""
},

{
    "location": "readings.html#",
    "page": "Readings",
    "title": "Readings",
    "category": "page",
    "text": "CurrentModule = SDDPOn this page, we've collated a variety of papers and books we think are helpful readings that cover knowledge needed to use SDDP.jl."
},

{
    "location": "readings.html#Stochastic-Optimization-1",
    "page": "Readings",
    "title": "Stochastic Optimization",
    "category": "section",
    "text": "A general primer on Stochastic ProgrammingBirge, J.R., Louveaux, F., 2011. Introduction to Stochastic Programming,  Springer Series in Operations Research and Financial Engineering. Springer New  York, New York, NY.  doi:10.1007/978-1-4614-0237-4Some overviews of Stochastic Optimization and where it sits in relation to other fieldsPowell, W.B., 2014. Clearing the Jungle of Stochastic Optimization, in: Newman,  A.M., Leung, J., Smith, J.C., Greenberg, H.J. (Eds.), Bridging Data and  Decisions. INFORMS, pp. 109–137. linkPowell, W.B., 2016. A Unified Framework for Optimization under Uncertainty  TutORials in Operations Research, in: Optimization Challenges in Complex,  Networked and Risky Systems. pp. 45–83. link"
},

{
    "location": "readings.html#Stochastic-Dual-Dynamic-Programming-1",
    "page": "Readings",
    "title": "Stochastic Dual Dynamic Programming",
    "category": "section",
    "text": "The original paper presenting SDDPPereira, M.V.F., Pinto, L.M.V.G., 1991. Multi-stage stochastic optimization  applied to energy planning. Mathematical Programming 52, 359–375. doi:10.1007/BF01582895The paper presenting the Markov version of SDDP implemented in this libraryPhilpott, A.B., de Matos, V.L., 2012. Dynamic sampling algorithms for multi-stage  stochastic programs with risk aversion. European Journal of Operational Research  218, 470–483. doi:10.1016/j.ejor.2011.10.056"
},

{
    "location": "readings.html#Julia-1",
    "page": "Readings",
    "title": "Julia",
    "category": "section",
    "text": "The organisation's websitehttps://julialang.org/The paper describing JuliaBezanson, J., Edelman, A., Karpinski, S., Shah, V.B., 2017. Julia: A Fresh  Approach to Numerical Computing. SIAM Review 59, 65–98. doi:10.1137/141000671"
},

{
    "location": "readings.html#JuMP-1",
    "page": "Readings",
    "title": "JuMP",
    "category": "section",
    "text": "Source code on Github[www.github.com/JuliaOpt/JuMP.jl]((https://www.github.com/JuliaOpt/JuMP.jl)The paper describing JuMPDunning, I., Huchette, J., Lubin, M., 2017. JuMP: A Modeling Language for  Mathematical Optimization. SIAM Review 59, 295–320. doi:10.1137/15M1020575"
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
    "location": "apireference.html#SDDP.@rhsnoise",
    "page": "Reference",
    "title": "SDDP.@rhsnoise",
    "category": "Macro",
    "text": "@rhsnoise(sp, rhs, constraint)\n\nDescription\n\nAdd a constraint with a noise in the RHS vector to the subproblem sp.\n\nArguments\n\nsp         the subproblem\nrhs        keyword argument key=value where value is a one-dimensional array containing the noise realisations\nconstraint any valid JuMP @constraint syntax that includes the keyword defined by rhs\n\nExamples\n\n@rhsnoise(sp, i=1:2, x + y <= i )\n@rhsnoise(sp, i=1:2, x + y <= 3 * rand(2)[i] )\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.@rhsnoises",
    "page": "Reference",
    "title": "SDDP.@rhsnoises",
    "category": "Macro",
    "text": "@rhsnoises(sp, rhs, begin\n    constraint\nend)\n\nDescription\n\nThe plural form of @rhsnoise similar to the JuMP macro @constraints.\n\nArguments\n\nSee @rhsnoise.\n\nExamples\n\n@rhsnoises(sp, i=1:2, begin\n    x + y <= i\n    x + y <= 3 * rand(2)[i]\nend)\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.setnoiseprobability!",
    "page": "Reference",
    "title": "SDDP.setnoiseprobability!",
    "category": "Function",
    "text": "setnoiseprobability!(sp::JuMP.Model, distribution::Vector{Float64})\n\nDescription\n\nSet the probability distribution of the stagewise independent noise in the sp subproblem.\n\nArguments\n\nsp            the subproblem\ndistribution vector containing the probability of each outcome occuring.   Should sum to 1. Defaults to the uniform distribution.\n\nExamples\n\nIf there are two realizations:\n\nsetnoiseprobability!(sp, [0.3, 0.7])\nsetnoiseprobability!(sp, [0.5, 0.6]) will error!\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.@stageobjective",
    "page": "Reference",
    "title": "SDDP.@stageobjective",
    "category": "Macro",
    "text": "@stageobjective!(sp, kw=noises, objective)\n\nDescription\n\nDefine an objective that depends on the realization of the stagewise noise. objective can be any valid third argument to the JuMP @objective macro (i.e. @objective(sp, Min, objective)) that utilises the variable kw that takes the realizations defined in noises.\n\nExamples\n\n@stageobjective(sp, w=1:2, w * x)\n@stageobjective(sp, i=1:2, w[i]^2 * x)\n@stageobjective(sp, i=1:2, x[i])\n\n\n\n"
},

{
    "location": "apireference.html#Communicating-the-problem-to-the-solver-1",
    "page": "Reference",
    "title": "Communicating the problem to the solver",
    "category": "section",
    "text": "SDDPModel\nExpectation\nNestedAVaR\nDefaultCutOracle\nLevelOneCutOracle\n@state\n@states\n@rhsnoise\n@rhsnoises\nsetnoiseprobability!\n@stageobjective"
},

{
    "location": "apireference.html#JuMP.solve",
    "page": "Reference",
    "title": "JuMP.solve",
    "category": "Function",
    "text": "solve(m::SDDPModel; kwargs...)\n\nDescription\n\nSolve the SDDPModel m using SDDP. Accepts a number of keyword arguments to control the solution process.\n\nPositional arguments\n\nm: the SDDPModel to solve\n\nKeyword arguments\n\nmax_iterations::Int:  The maximum number of cuts to add to a single stage problem before terminating.  Defaults to 10.\ntime_limit::Real:  The maximum number of seconds (in real time) to compute for before termination.  Defaults to Inf.\nsimulation::MonteCarloSimulation: see MonteCarloSimulation\nbound_convergence::BoundConvergence: see BoundConvergence\ncut_selection_frequency::Int:  Frequency (by iteration) with which to rebuild subproblems using a subset of  cuts. Frequent cut selection (i.e. cut_selection_frequency is small) reduces  the size of the subproblems that are solved, but incurrs the overhead of rebuilding  the subproblems. However, infrequent cut selection (i.e.  cut_selection_frequency is large) allows the subproblems to grow large (many  constraints) leading to an increase in the solve time of individual subproblems.  Defaults to 0 (never run).\nprint_level::Int:   0 - off: nothing logged.   1 - on: solve iterations logged.   2 - verbose: detailed timing information is also logged.   Defaults to 1\nlog_file::String:  Relative filename to write the log to disk. Defaults to \"\" (no log written)\nsolve_type:  One of\nAsyncronous() - solve using a parallelised algorithm\nSerial() - solve using a serial algorithm\nDefault chooses automatically based on the number of available processors.\nreduce_memory_footprint::Bool:  Implements the idea proposed in https://github.com/JuliaOpt/JuMP.jl/issues/969#issuecomment-282191105  to reduce the memory consumption when running SDDP. This is an issue if you  wish to save the model m to disk since it discards important information.  Defaults to false.\ncut_output_file::String:  Relative filename to write discovered cuts to disk. Defaults to \"\" (no cuts written)\n\nReturns\n\nstatus::Symbol:  Reason for termination. One of\n:solving\n:interrupted\n:converged\n:max_iterations\n:bound_convergence\n:time_limit\n\n\n\n"
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
    "location": "apireference.html#Solving-the-problem-efficiently-1",
    "page": "Reference",
    "title": "Solving the problem efficiently",
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
    "text": " SDDP.plotvaluefunction(m::SDDPModel, stage::Int, markovstate::Int, states::Union{Float64, AbstractVector{Float64}}...; label1=\"State 1\", label2=\"State 2\")\n\nDescription\n\nPlot the value function of stage stage and Markov state markovstate in the SDDPModel m at the points in the discretized state space given by states. If the value in states is a real number, the state is evaluated at that point. If the value is a vector, the state is evaluated at all the points in the vector. At most two states can be vectors.\n\nExamples\n\nSDDP.plotvaluefunction(m, 2, 1, 0.0:0.1:1.0, 0.5, 0.0:0.1:1.0; label1=\"State 1\", label2=\"State 3\")\n\n\n\n"
},

{
    "location": "apireference.html#Understanding-the-solution-1",
    "page": "Reference",
    "title": "Understanding the solution",
    "category": "section",
    "text": "simulate\ngetbound\nnewplot\naddplot!\nshow\nplotvaluefunction"
},

{
    "location": "apireference.html#SDDP.loadcuts!",
    "page": "Reference",
    "title": "SDDP.loadcuts!",
    "category": "Function",
    "text": "loadcuts!(m::SDDPModel, filename::String)\n\nLoad cuts from the file created using the cut_output_file argument in solve.\n\nExample\n\nm = SDDPModel() do ... end\nstatus = solve(m; cut_output_file=\"path/to/m.cuts\")`\nm2 = SDDPModel() do ... end\nloadcuts!(m2, \"path/to/m.cuts\")\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.savemodel!",
    "page": "Reference",
    "title": "SDDP.savemodel!",
    "category": "Function",
    "text": "SDDP.savemodel!(filename::String, m::SDDPModel)\n\nSave the SDDPModel m to the location filename. Can be loaded at a later date with m = SDDP.loadmodel(filename).\n\nNote: this function relies in the internal Julia Base.serializefunction. It should not be relied on to save an load models between versions of Julia (i.e between v0.5 and v0.6). For a longer-term solution, see SDDP.loadcuts! for help.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.loadmodel",
    "page": "Reference",
    "title": "SDDP.loadmodel",
    "category": "Function",
    "text": "loadmodel(filename::String)\n\nLoad a model from the location filename that was saved using SDDP.savemodel!.\n\nNote: this function relies in the internal Julia Base.serializefunction. It should not be relied on to save an load models between versions of Julia (i.e between v0.5 and v0.6). For a longer-term solution, see SDDP.loadcuts! for help.\n\n\n\n"
},

{
    "location": "apireference.html#Read-and-write-the-model-to-disk-1",
    "page": "Reference",
    "title": "Read and write the model to disk",
    "category": "section",
    "text": "    loadcuts!\n    savemodel!\n    loadmodel"
},

{
    "location": "apireference.html#SDDP.AbstractCutOracle",
    "page": "Reference",
    "title": "SDDP.AbstractCutOracle",
    "category": "Type",
    "text": "AbstractCutOracle\n\nDescription\n\nAbstract type for all cut oracles.\n\n\n\n"
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
    "location": "apireference.html#SDDP.AbstractRiskMeasure",
    "page": "Reference",
    "title": "SDDP.AbstractRiskMeasure",
    "category": "Type",
    "text": "AbstractRiskMeasure\n\nDescription\n\nAbstract type for all risk measures.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.modifyprobability!",
    "page": "Reference",
    "title": "SDDP.modifyprobability!",
    "category": "Function",
    "text": "modifyprobability!(measure::AbstractRiskMeasure,\n        riskadjusted_distribution,\n        original_distribution::Vector{Float64},\n        observations::Vector{Float64},\n        m::SDDPModel,\n        sp::JuMP.Model\n)\n\nDescription\n\nCalculate the risk-adjusted probability of each scenario using the 'change-of-probabilities' approach of Philpott, de Matos, and Finardi,(2013). On solving multistage stochastic programs with coherent risk measures. Operations Research 61(4), 957-970.\n\nArguments\n\nmeasure::AbstractRiskMeasure\n\nThe risk measure\n\nriskadjusted_distribution\n\nA new probability distribution\n\noriginal_distribution::Vector{Float64}\n\nThe original probability distribution.\n\nobservations::Vector{Float64}\n\nThe vector of objective values from the next stage  problems (one for each scenario).\n\nm::SDDPModel\n\nThe full SDDP model\n\nsp::JuMP.Model\n\nThe stage problem that the cut will be added to.\n\n\n\n"
},

{
    "location": "apireference.html#Extras-for-Experts-1",
    "page": "Reference",
    "title": "Extras for Experts",
    "category": "section",
    "text": "AbstractCutOracle\nstorecut!\nvalidcutsAbstractRiskMeasure\nmodifyprobability!"
},

]}
