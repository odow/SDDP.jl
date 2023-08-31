var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "CurrentModule = SDDP"
},

{
    "location": "index.html#SDDP.jl-Documentation-1",
    "page": "Home",
    "title": "SDDP.jl Documentation",
    "category": "section",
    "text": "SDDP.jl is a package for solving large multistage convex stochastic optimization problems using stochastic dual dynamic programming. In this manual, we\'re going to assume a reasonable amount of background knowledge about stochastic optimization, the SDDP algorithm, Julia, and JuMP.note: Note\nIf you don\'t have that background, you may want to brush up on some Readings. Part I of my thesis may also be useful as it contains a primer on how to formulate multistage stochastic optimization problems (Chapter One), as well as an introduction and literature review of the SDDP algorithm (Chapter Two)."
},

{
    "location": "index.html#Getting-started-1",
    "page": "Home",
    "title": "Getting started",
    "category": "section",
    "text": "This package is unregistered so you will need to Pkg.clone it as follows:Pkg.clone(\"https://github.com/odow/SDDP.jl.git\")If you want to use the parallel features of SDDP.jl, you should start Julia with some worker processes (julia -p N), or add by running julia> addprocs(N) in a running Julia session.Once you\'ve got SDDP.jl installed, you should read some tutorials, beginning with Tutorial One: first steps."
},

{
    "location": "index.html#Citing-SDDP.jl-1",
    "page": "Home",
    "title": "Citing SDDP.jl",
    "category": "section",
    "text": "If you use SDDP.jl, we ask that you please cite the following paper:@article{dowson_sddp.jl,\n	title = {{SDDP}.jl: a {Julia} package for stochastic dual dynamic programming},\n	url = {http://www.optimization-online.org/DB_HTML/2017/12/6388.html},\n	journal = {Optimization Online},\n	author = {Dowson, Oscar and Kapelevich, Lea},\n	year = {2017}\n}"
},

{
    "location": "index.html#FAQ-1",
    "page": "Home",
    "title": "FAQ",
    "category": "section",
    "text": "Q. How do I make the constraint coefficients random?A. Due to the design of JuMP, it\'s difficult to efficiently modify constraint coefficients. Therefore, you can only vary the right hand-side of a constraint using the @rhsnoise macro.As a work around, we suggest you either reformulate the model so the uncertainty appears in the RHS, or model the uncertainty as a Markov process. Tutorial Four: Markovian policy graphs explains how to implement this. You might also want to take a look at the asset management example to see an example of this. Make sure you keep in mind that a new value function is built at each Markov state which increases the computation time and memory requirements."
},

{
    "location": "tutorial/01_first_steps.html#",
    "page": "Tutorial One: first steps",
    "title": "Tutorial One: first steps",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/01_first_steps.html#Tutorial-One:-first-steps-1",
    "page": "Tutorial One: first steps",
    "title": "Tutorial One: first steps",
    "category": "section",
    "text": "Hydrothermal scheduling is the most common application of stochastic dual dynamic programming. To illustrate some of the basic functionality of SDDP.jl, we implement a very simple model of the hydrothermal scheduling problem.In this model, there are two generators: a thermal generator, and a hydro generator. The thermal generator has a short-run marginal cost of \\$50/MWh in the first stage, \\$100/MWh in the second stage, and \\$150/MWh in the third stage. The hydro generator has a short-run marginal cost of \\$0/MWh.We consider the problem of scheduling the generation over three time periods in order to meet a known demand of 150 MWh in each period.The hydro generator draws water from a reservoir which has a maximum capacity of 200 units. We assume that at the start of the first time period, the reservoir is full. In addition to the ability to generate electricity by passing water through the hydroelectric turbine, the hydro generator can also spill water down a spillway (bypassing the turbine) in order to prevent the water from over-topping the dam. We assume that there is no cost of spillage.The objective of the optimization is to minimize the expected cost of generation over the three time periods."
},

{
    "location": "tutorial/01_first_steps.html#Formulating-the-problem-1",
    "page": "Tutorial One: first steps",
    "title": "Formulating the problem",
    "category": "section",
    "text": "First, we need to load some packages. For this example, we are going to use the Clp.jl package; however, you are free to use any solver that you could normally use with JuMP.using SDDP, JuMP, ClpNext, we need to initialize our model. In our example, we are minimizing, there are three stages, and we know a lower bound of 0.0. Therefore, we can initialize our model using the SDDPModel constructor:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0\n                                        ) do sp, t\n    # ... stuff to go here ...\nendIf you haven\'t seen the do sp, t ... end syntax before, this syntax is equivalent to the following:function define_subproblem(sp::JuMP.Model, t::Int)\n    # ... stuff to go here ...\nend\nm = SDDPModel(\n    define_subproblem,\n    sense           = :Min,\n    stages          = 3,\n    solver          = ClpSolver(),\n    objective_bound = 0.0\n)The function define_subproblem (although you can call it anything you like) takes two arguments: sp, a JuMP.Model that we will use to build each subproblem; and t, an Int that is a counter from 1 to the number of stages. In this case, t=1, 2, 3. The sense, stages, and solver keyword arguments to SDDPModel should be obvious; however, the objective_bound is worth explaining.In order to solve a model using SDDP, we need to define a valid lower bound for every subproblem. (See Introduction to SDDP for details.) In this example, the least-cost solution is to meet demand entirely from the hydro generator, incurring a cost of \\$0/MWh. Therefore, we set objective_bound=0.0.Now we need to define build each subproblem using a mix of JuMP and SDDP.jl syntax."
},

{
    "location": "tutorial/01_first_steps.html#State-variables-1",
    "page": "Tutorial One: first steps",
    "title": "State variables",
    "category": "section",
    "text": "There is one state variable in our model: the quantity of water in the reservoir at the end of stage t. Two add this state variable to the model, SDDP.jl defines the @state macro.  This macro takes three arguments:sp - the JuMP model;\nan expression for the outgoing state variable; and\nan expression for the incoming state variable.The 2nd argument can be any valid JuMP @variable syntax and can include, for example, upper and lower bounds. The 3rd argument must be the name of the incoming state variable, followed by ==, and then the value of the state variable at the root node of the policy graph. For our hydrothermal example, the state variable can be constructed as:@state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)note: Note\nYou must define the same state variables (in the same order) in every subproblem. If a state variable is unused in a particular stage, it must still be defined."
},

{
    "location": "tutorial/01_first_steps.html#Control-variables-1",
    "page": "Tutorial One: first steps",
    "title": "Control variables",
    "category": "section",
    "text": "We now need to define some control variables. In SDDP.jl, control variables are just normal JuMP variables. Therefore, we can define the three variables in the hydrothermal scheduling problem (thermal generation, hydro generation, and the quantity of water to spill) as follows:@variables(sp, begin\n    thermal_generation >= 0\n    hydro_generation   >= 0\n    hydro_spill        >= 0\n end)"
},

{
    "location": "tutorial/01_first_steps.html#Constraints-1",
    "page": "Tutorial One: first steps",
    "title": "Constraints",
    "category": "section",
    "text": "Before we specify the constraints, we need to create some data. For this problem, we need the inflow to the reservoir in each stage t=1, 2, 3. Therefore, we create the vector:inflow = [50.0, 50.0, 50.0]The inflow in stage t can be accessed as inflow[t].First, we have the water balance constraint: the volume of water at the end of the stage must equal the volume of water at the start of the stage, plus any inflows, less that used for generation or spilled down the spillway.@constraint(sp,\n    incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume\n)Note that we use t defined by the SDDPModel constructor. There is also a constraint that total generation must equal demand of 150 MWh:@constraint(sp,\n    thermal_generation + hydro_generation == 150\n)"
},

{
    "location": "tutorial/01_first_steps.html#The-stage-objective-1",
    "page": "Tutorial One: first steps",
    "title": "The stage objective",
    "category": "section",
    "text": "Finally, there is a cost on thermal generation of \\$50/MWh in the first stage, \\$100/MWh in the second stage, and \\$150/MWh in the third stage. To add the stage-objective, we use the aptly named @stageobjective macro provided by SDDP.jl:if t == 1\n    @stageobjective(sp,  50.0 * thermal_generation )\nelseif t == 2\n    @stageobjective(sp, 100.0 * thermal_generation )\nelseif t == 3\n    @stageobjective(sp, 150.0 * thermal_generation )\nendinfo: Info\nif statements can be used more broadly in the subproblem definition to conditionally and variables and constraints into different subproblems.We can also implement the stage-objective more succinctly using a vector:fuel_cost = [50.0, 100.0, 150.0]\n@stageobjective(sp, fuel_cost[t] * thermal_generation )"
},

{
    "location": "tutorial/01_first_steps.html#Solving-the-problem-1",
    "page": "Tutorial One: first steps",
    "title": "Solving the problem",
    "category": "section",
    "text": "Putting all that we have discussed above together, we get:using SDDP, JuMP, Clp\nm = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0\n                                        ) do sp, t\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n     end)\n    inflow = [50.0, 50.0, 50.0]\n    @constraints(sp, begin\n        incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, fuel_cost[t] * thermal_generation )\nendTo solve this problem, we use the solve method:status = solve(m; iteration_limit=5)The argument iteration_limit is self-explanatory. The return value status is a symbol describing why the SDDP algorithm terminated. In this case, the value is :iteration_limit. We discuss other arguments to the solve method and other possible values for status in future sections of this manual.During the solve, the following log is printed to the screen.-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Serial solver\n    Model:\n        Stages:         3\n        States:         1\n        Subproblems:    3\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n       15.000K         5.000K        |     1    0.0      0    0.0    0.0\n        5.000K         5.000K        |     2    0.0      0    0.0    0.0\n        5.000K         5.000K        |     3    0.0      0    0.0    0.0\n        5.000K         5.000K        |     4    0.0      0    0.0    0.0\n        5.000K         5.000K        |     5    0.0      0    0.0    0.0\n-------------------------------------------------------------------------------\n    Other Statistics:\n        Iterations:         5\n        Termination Status: iteration_limit\n===============================================================================The header and footer of the output log contain self-explanatory statistics about the problem. The numeric columns are worthy of description. Each row corresponds to one iteration of the SDDP algorithm.The left half of the log relates to the objective of the problem. In the Simulation column, we give the cumulative cost of each forward pass. In the Bound column, we give the lower bound (upper if maximizing) obtained after the backward pass has completed in each iteration. Ignore the % Gap column for now, that is addressed in Tutorial Two: RHS noise.The right half of the log displays timing statistics. Cut Passes displays the number of cutting iterations conducted (in #) and the time it took to (in Time). Ignore the Simulations columns for now, they are addressed in Tutorial Tutorial Two: RHS noise. Finally, the Total Time column records the total time spent solving the problem.This log can be silenced by setting the print_level keyword argument to solve to 0. In addition, the log will be written to the file given by the log_file keyword argument (this is off by default)."
},

{
    "location": "tutorial/01_first_steps.html#Understanding-the-solution-1",
    "page": "Tutorial One: first steps",
    "title": "Understanding the solution",
    "category": "section",
    "text": "The first thing we want to do is to query the lower (upper if maximizing) bound of the solution. This can be done via the getbound function:getbound(m)This returns the value of the Bound column in the last row in the output table above. In this example, the bound is 5000.0.Then, we can perform a Monte Carlo simulation of the policy using the simulate function. It takes three arguments. The first is the SDDPModel m. The second is the number of replications to perform. The third is a vector of variable names to record the value of at each stage and replication. Since our example is deterministic, it is sufficient to perform a single replication:simulation_result = simulate(m,\n    1,\n    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill]\n)The return value, simulation_result, is a vector of dictionaries containing one element for each Monte Carlo replication. In this case, length(simulation_result) = 1. The keys of the dictionary are the variable symbols given in the simulate function, and their associated values are vectors, with one element for each stage, or the variable value in the simulated solution. For example, we can query the optimal quantity of hydro generation in each stage as follows:julia> simulation_result[1][:hydro_generation]\n3-element Array{Any, 1}:\n  50.0\n 150.0\n 150.0This concludes our first very simple tutorial for SDDP.jl. In the next tutorial, Tutorial Two: RHS noise, we introduce stagewise-independent noise into the model."
},

{
    "location": "tutorial/02_rhs_noise.html#",
    "page": "Tutorial Two: RHS noise",
    "title": "Tutorial Two: RHS noise",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/02_rhs_noise.html#Tutorial-Two:-RHS-noise-1",
    "page": "Tutorial Two: RHS noise",
    "title": "Tutorial Two: RHS noise",
    "category": "section",
    "text": "In Tutorial One: first steps, we formulated a simple hydrothermal scheduling problem. In this tutorial, we extend the model to include stagewise-independent noise in the right-hand side of the constraints.note: Note\nNotably, SDDP.jl does not allow stagewise-independent noise terms in the constraint matrix. However, this can be modelled using a Markovian policy graph like the one in Tutorial Four: Markovian policy graphs.Recall that our model for the hydrothermal scheduling problem  from Tutorial One: first steps is:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0\n                                        ) do sp, t\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n     end)\n    inflow = [50.0, 50.0, 50.0]\n    @constraints(sp, begin\n        incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, fuel_cost[t] * thermal_generation )\nend"
},

{
    "location": "tutorial/02_rhs_noise.html#Formulating-the-problem-1",
    "page": "Tutorial Two: RHS noise",
    "title": "Formulating the problem",
    "category": "section",
    "text": "In this tutorial, we are going to model inflows that are stagewise-independent. Specifically, we assume that in each stage, there is an even probability of sampling an inflow of 0.0, 50.0, or 100.0. To add this noise term to the model, we need to use the @rhsnoise macro provided by SDDP.jl.@rhsnoise is similar to the JuMP @constraint macro. It takes three arguments. The first is the subproblem sp. The second argument is of the form name=[realizations], where name is a descriptive name, and realizations is a vector of elements in the sample space. The third argument is any valid JuMP constraint that utilizes name in the right-hand side. For our example, we have:@rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n    outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n)However, the realizations do not have to be the full right-hand side term. The following is also valid:inflows = [0.0, 50.0, 100.0]\n@rhsnoise(sp, i = [1,2,3],\n    outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflows[i]\n)We can set the probability of sampling each element in the sample space using the setnoiseprobability! function. If setnoiseprobability! isn\'t called, the distribution is assumed to be uniform. Despite this, for the sake of completeness, we set the probability for our example as:setnoiseprobability!(sp, [1/3, 1/3, 1/3])Our model is now:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0\n                                        ) do sp, t\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    setnoiseprobability!(sp, [1/3, 1/3, 1/3])\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, fuel_cost[t] * thermal_generation )\nend"
},

{
    "location": "tutorial/02_rhs_noise.html#Solving-the-problem-1",
    "page": "Tutorial Two: RHS noise",
    "title": "Solving the problem",
    "category": "section",
    "text": "Now we need to solve the problem. As in Tutorial One: first steps, we use the solve function. However, this time we utilize some additional arguments.Since our problem is stochastic, we often want to simulate the policy in order to estimate the upper (lower if maximizing) bound. This can be controlled via the simulation keyword to solve. The syntax has a lot going on so we\'re going to give an example of how it is used, and then walk through the different components.status = solve(m,\n    simulation = MonteCarloSimulation(\n        frequency  = 2,\n        confidence = 0.95,\n        terminate  = true\n        min        = 50,\n        step       = 50,\n        max        = 100,\n    )\n)First, the frequency argument specifies how often the  Monte Carlo simulation is conducted (iterations/simulation). For this example, we conduct a Monte Carlo simulation every two iterations. Second, the confidence specifies the level at which to conduct the confidence interval. In this example, we construct a 95% confidence interval. Third, the terminate argument is a Boolean defining if we should terminate the method if the lower limit of the confidence interval is less than the lower bound (upper limit and bound for maximization problems). The final three arguments implement the method of sequential sampling: min gives the minimum number of replications to conduct before the construction of a confidence interval. If there is evidence of convergence, another step replications are conducted. This continues until either: (1) max number of replications have been conducted; or (2) there is no evidence of convergence. For our example, we conduct 50 replications, and if there is no evidence of convergence, we conduct another 50 replications.The output from the log is now:-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Serial solver\n    Model:\n        Stages:         3\n        States:         1\n        Subproblems:    3\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n       17.500K         3.438K        |     1    0.0      0    0.0    0.0\n   7.606K   10.894K    7.500K   1.4  |     2    0.0     50    0.0    0.0\n        7.500K         8.333K        |     3    0.0     50    0.0    0.0\n   7.399K    9.651K    8.333K -11.2  |     4    0.0    150    0.1    0.1\n-------------------------------------------------------------------------------\n    Other Statistics:\n        Iterations:         4\n        Termination Status: converged\n===============================================================================Compared with the log of a solve without using the simulation keyword, a few things have changed. First, in the second and fourth rows (i.e. the iterations in which a Monte Carlo simulation was conducted) the Simulation column now gives two values. This pair is the confidence interval for the estimate of the upper bound.Second, in iterations in which a Monte Carlo simulation is conducted, there is an entry in the % Gap column. This gaps measures the difference between the lower limit of the simulated confidence interval and the lower bound (given in the Bound) column. If the gap is positive, there is evidence that the model has not converged. Once the gap is negative, the lower bound lies above the lower limit of the confidence interval and we can terminate the algorithm.The third difference is that the Simulations column now records the number of Monte Replications conducted to estimate the upper bound (in #) and time performing those Monte Carlo replications (in Time). You can use this information to tune the frequency at which the policy is tested for convergence.Also observe that the first time we performed the Monte Carlo simulation, we only conducted 50 replications; however, the second time we conducted 100. This demonstrates the sequential sampling method at work.Finally, the termination status is now :converged instead of :iteration_limit."
},

{
    "location": "tutorial/02_rhs_noise.html#Understanding-the-solution-1",
    "page": "Tutorial Two: RHS noise",
    "title": "Understanding the solution",
    "category": "section",
    "text": "The first thing we want to do is to query the lower (upper if maximizing) bound of the solution. This can be done via the getbound function:getbound(m)This returns the value of the Bound column in the last row in the output table above. In this example, the bound is 8333.0.Then, we can perform a Monte Carlo simulation of the policy using the simulate function. We perform 500 replications and record the same variables as we did in Tutorial One: first steps.simulation_result = simulate(m,\n    500,\n    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill]\n)This time, length(simulation_result) = 500. In addition to the variables, we also record some additional fields. This includes :stageobjective, the value of the stage-objective in each stage. We can calculate the cumulative objective of each replication by summing the stage-objectives as follows:julia> sum(simulation_result[100][:stageobjective])\n2500.0We can calculate the objective of all of each replication using Julia\'s generator syntax:julia> objectives = [sum(replication[:stageobjective]) for replication in simulation_result]\n500-element Array{Float64, 1}:\n  5000.0\n 20000.0\n 15000.0\n ⋮Then, we can calculate the mean and standard deviation of these objectives:julia> mean(objectives), std(objectives)\n(8025.0, 5567.66)We can query the noise that was sampled in each stage using the :noise key. This returns the index of the noise from the vector of realizations. For example:julia> simulation_result[100][:noise]\n3-element Array{Int, 1}:\n 1\n 3\n 2This concludes our second tutorial for SDDP.jl. In the next tutorial, Tutorial Three: objective noise, we introduce stagewise-independent noise into the objective function."
},

{
    "location": "tutorial/03_objective_noise.html#",
    "page": "Tutorial Three: objective noise",
    "title": "Tutorial Three: objective noise",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/03_objective_noise.html#Tutorial-Three:-objective-noise-1",
    "page": "Tutorial Three: objective noise",
    "title": "Tutorial Three: objective noise",
    "category": "section",
    "text": "In Tutorial One: first steps, we formulated a simple, deterministic hydrothermal scheduling problem. Then, Tutorial Two: RHS noise, we extended this model to include stagewise-independent noise to the right-hand side of a constraint. Now, in this tutorial, we extend the model to include stagewise-independent noise in the objective function.note: Note\nNotably, SDDP.jl does not allow stagewise-independent noise terms in the constraint matrix. However, this can be modelled using a Markovian policy graph like the one in Tutorial Four: Markovian policy graphs.Recall that our model for the hydrothermal scheduling problem  from Tutorial Two: RHS noise is:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0\n                                        ) do sp, t\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    setnoiseprobability!(sp, [1/3, 1/3, 1/3])\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, fuel_cost[t] * thermal_generation )\nend"
},

{
    "location": "tutorial/03_objective_noise.html#Formulating-the-problem-1",
    "page": "Tutorial Three: objective noise",
    "title": "Formulating the problem",
    "category": "section",
    "text": "In this tutorial, we are going to model the fuel cost of the thermal generator by a stagewise-independent process. Specifically, we assume that in each stage, there is an even probability of sampling a fuel cost of 80%, 100%, or 120% of the usual fuel costs of \\$50/MWh in the first stage, \\$100/MWh in the second stage, and \\$150/MWh in the third stage. To add this noise term to the model, we need to use a modified version of the @stageobjective macro provided by SDDP.jl.This version of @stageobjective is similar to the @rhsnoise macro that we discussed in Tutorial Two: RHS noise. It takes three arguments. The first is the subproblem sp. The second argument is of the form name=[realizations], where name is a descriptive name, and realizations is a vector of elements in the sample space. The third argument is any valid input to the normal @stageobjective macro.It is important to note that there must be the same number of realizations in the objective as there are realizations in the right-hand-side random variable (created using @rhsnoise). The two noise terms will be sampled in unison, so that when the first element of the right-hand side noise is sampled, so to will the first element of the objective noise. If the two noise terms should be sampled independently, the user should form the Cartesian product.For our example, the price multiplier and the inflows are negatively correlated. Therefore, when inflow=0.0, the multiplier is 1.2, when inflow=50.0, the multiplier is 1.0, and when inflow=100.0, the multiplier is 0.8. Thus, we have:fuel_cost = [50.0, 100.0, 150.0]\n@stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],\n    mupliplier * fuel_cost[t] * thermal_generation\n)As in Tutorial Two: RHS noise, the noise terms are sampled using the probability distribution set by the setnoiseprobability! function.Our model is now:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0\n                                        ) do sp, t\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    setnoiseprobability!(sp, [1/3, 1/3, 1/3])\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],\n        mupliplier * fuel_cost[t] * thermal_generation\n    )\nend"
},

{
    "location": "tutorial/03_objective_noise.html#Solving-the-problem-1",
    "page": "Tutorial Three: objective noise",
    "title": "Solving the problem",
    "category": "section",
    "text": "Now we need to solve the problem. As in the previous two tutorials, we use the solve function. However, this time we use the bound stalling stopping rule. This can be controlled via the bound_stalling keyword to solve. The syntax has a lot going on so we\'re going to give an example of how it is used, and then walk through the different components.status = solve(m,\n    bound_stalling = BoundStalling(\n        iterations = 5,\n        rtol       = 0.0,\n        atol       = 1e-6\n    )\n)First, the iterations argument specifies how many iterations that the bound must change by less than atol or rtol before terminating. For this example, we choose to terminate the SDDP algorithm after the bound has failed to improve for 5 iterations. Second, the rtol and atol keywords determine the absolute and relative tolerances by which we compare the equality of the bound in consecutive iterations. In this case, since the model is simple we choose an absolute convergence tolerance of 1e-6.The termination status is :bound_stalling, and the output from the log is now:-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Serial solver\n    Model:\n        Stages:         3\n        States:         1\n        Subproblems:    3\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n        7.500K         5.733K        |     1    0.0      0    0.0    0.0\n       11.800K         8.889K        |     2    0.0      0    0.0    0.0\n       14.000K         9.167K        |     3    0.0      0    0.0    0.0\n       11.000K         9.167K        |     4    0.0      0    0.0    0.0\n        9.000K         9.167K        |     5    0.0      0    0.0    0.0\n        2.000K         9.167K        |     6    0.0      0    0.0    0.0\n       14.000K         9.167K        |     7    0.0      0    0.0    0.0\n-------------------------------------------------------------------------------\n    Other Statistics:\n        Iterations:         7\n        Termination Status: bound_stalling\n==============================================================================="
},

{
    "location": "tutorial/03_objective_noise.html#Understanding-the-solution-1",
    "page": "Tutorial Three: objective noise",
    "title": "Understanding the solution",
    "category": "section",
    "text": "Instead of performing a Monte Carlo simulation, you may want to simulate one particular sequence of noise realizations. This historical simulation can also be conducted using the simulate function.simulation_result = simulate(m,\n    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill],\n    noises = [1, 1, 3]\n)This time, simulation_result is a single dictionary. We can query the objective of the simulation as follows:julia> simulation_result[:objective]\n9000.0Interestingly, despite sampling the low-inflow, high-price realization in the first stage, the model generates 150 MWh at a price of \\$60/MWh:julia> simulation_result[:thermal_generation]\n3-element Array{Any, 1}:\n 150.0\n   0.0\n   0.0This concludes our third tutorial for SDDP.jl. In the next tutorial, Tutorial Four: Markovian policy graphs, we introduce stagewise-dependent noise via a Markov chain."
},

{
    "location": "tutorial/04_markovian_policygraphs.html#",
    "page": "Tutorial Four: Markovian policy graphs",
    "title": "Tutorial Four: Markovian policy graphs",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/04_markovian_policygraphs.html#Tutorial-Four:-Markovian-policy-graphs-1",
    "page": "Tutorial Four: Markovian policy graphs",
    "title": "Tutorial Four: Markovian policy graphs",
    "category": "section",
    "text": "In our three tutorials (Tutorial One: first steps, Tutorial Two: RHS noise, and Tutorial Three: objective noise), we formulated a simple hydrothermal scheduling problem with stagewise-independent noise in the right-hand side of the constraints and in the objective function. Now, in this tutorial, we introduce some stagewise-dependent uncertainty using a Markov chain.Recall that our model for the hydrothermal scheduling problem  from Tutorial Three: objective noise is:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0\n                                        ) do sp, t\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    setnoiseprobability!(sp, [1/3, 1/3, 1/3])\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],\n        mupliplier * fuel_cost[t] * thermal_generation\n    )\nend"
},

{
    "location": "tutorial/04_markovian_policygraphs.html#Formulating-the-problem-1",
    "page": "Tutorial Four: Markovian policy graphs",
    "title": "Formulating the problem",
    "category": "section",
    "text": "In this tutorial we consider a Markov chain with two climate states: wet and dry. Each Markov state is associated with an integer, in this case the wet climate state  is Markov state 1 and the dry climate state is Markov state 2. In the wet climate state, the probability of the high inflow increases to 50%, and the probability of the low inflow decreases to 1/6. In the dry climate state, the converse happens. There is also persistence in the climate state: the probability of remaining in the current state is 75%, and the probability of transitioning to the other climate state is 25%. We assume that the first stage starts in the wet climate state.For each stage, we need to provide a Markov transition matrix. This is an MxN matrix, where the element A[i,j] gives the probability of transitioning from Markov state i in the previous stage to Markov state j in the current stage. The first stage is special because we assume there is a \"zero\'th\" stage which has one Markov state. Furthermore, the number of columns in the transition matrix of a stage (i.e. the number of Markov states) must equal the number of rows in the next stage\'s transition matrix. For our example, the vector of Markov transition matrices is given by:T = Array{Float64, 2}[\n    [ 1.0 0.0 ],\n    [ 0.75 0.25 ; 0.25 0.75 ],\n    [ 0.75 0.25 ; 0.25 0.75 ]\n]However, note that we never sample the dry Markov state in stage one. Therefore, we can drop that Markov state so that there is only one Markov state in stage 1. We also need to modify the transition matrix in stage 2 to account for this:T = Array{Float64, 2}[\n    [ 1.0 ]\',\n    [ 0.75 0.25 ],\n    [ 0.75 0.25 ; 0.25 0.75 ]\n]To add the Markov chain to the model, we modifications are required. First, we give the vector of transition matrices to the SDDPModel constructor using the markov_transition keyword. Second, the do sp, t ... end syntax is extended to do sp, t, i ... end, where i is the index of the Markov state and runs from i=1 to the number of Markov states in stage t. Now, both t and i can be used anywhere inside the subproblem definition.Our model is now:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0,\n      markov_transition = Array{Float64, 2}[\n          [ 1.0 ]\',\n          [ 0.75 0.25 ],\n          [ 0.75 0.25 ; 0.25 0.75 ]\n      ]\n                                        ) do sp, t, i\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],\n        mupliplier * fuel_cost[t] * thermal_generation\n    )\n    if i == 1  # wet climate state\n        setnoiseprobability!(sp, [1/6, 1/3, 0.5])\n    else       # dry climate state\n        setnoiseprobability!(sp, [0.5, 1/3, 1/6])\n    end\nend"
},

{
    "location": "tutorial/04_markovian_policygraphs.html#Solving-the-problem-1",
    "page": "Tutorial Four: Markovian policy graphs",
    "title": "Solving the problem",
    "category": "section",
    "text": "Now we need to solve the problem. As in the previous two tutorials, we use the solve function. However, this time we terminate the SDDP algorithm by setting a time limit (in seconds) using the time_limit keyword:status = solve(m,\n    time_limit = 0.05\n)The termination status is :time_limit, and the output from the log is now:-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Serial solver\n    Model:\n        Stages:         3\n        States:         1\n        Subproblems:    3\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n        0.000          6.198K        |     1    0.0      0    0.0    0.0\n        2.000K         7.050K        |     2    0.0      0    0.0    0.0\n        2.000K         7.050K        |     3    0.0      0    0.0    0.0\n        2.000K         7.135K        |     4    0.0      0    0.0    0.0\n        5.000K         7.135K        |     5    0.0      0    0.0    0.0\n        2.000K         7.135K        |     6    0.0      0    0.0    0.0\n        2.000K         7.135K        |     7    0.0      0    0.0    0.0\n        2.000K         7.135K        |     8    0.0      0    0.0    0.0\n        2.000K         7.135K        |     9    0.0      0    0.0    0.0\n        9.000K         7.135K        |    10    0.0      0    0.0    0.0\n        2.000K         7.135K        |    11    0.0      0    0.0    0.0\n        5.000K         7.135K        |    12    0.1      0    0.0    0.1\n-------------------------------------------------------------------------------\n    Other Statistics:\n        Iterations:         12\n        Termination Status: time_limit\n==============================================================================="
},

{
    "location": "tutorial/04_markovian_policygraphs.html#Understanding-the-solution-1",
    "page": "Tutorial Four: Markovian policy graphs",
    "title": "Understanding the solution",
    "category": "section",
    "text": "Instead of performing a Monte Carlo simulation, you may want to simulate one particular sequence of noise realizations. This historical simulation can also be conducted using the simulate function.simulation_result = simulate(m,\n    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill],\n    noises       = [1, 1, 3],\n    markovstates = [1, 2, 2]\n)Again, simulation_result is a single dictionary. In addition to the variable values and the special keys :noise, :objective, and :stageobjective, SDDP.jl also records the index of the Markov state in each stage via the :markov key. We can confirm that the historical sequence of Markov states was visited as follows:julia> simulation_result[:markov]\n3-element Array{Int, 1}:\n 1\n 2\n 2This concludes our fourth tutorial for SDDP.jl. In the next tutorial, Tutorial Five: risk, we introduce risk into the problem."
},

{
    "location": "tutorial/05_risk.html#",
    "page": "Tutorial Five: risk",
    "title": "Tutorial Five: risk",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/05_risk.html#Tutorial-Five:-risk-1",
    "page": "Tutorial Five: risk",
    "title": "Tutorial Five: risk",
    "category": "section",
    "text": "Over the previous four tutorials, we formulated a simple hydrothermal scheduling problem. Now, in this tutorial, we introduce some risk into the model using nested risk measures.Recall that our model for the hydrothermal scheduling problem from Tutorial Four: Markovian policy graphs is:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0,\n      markov_transition = Array{Float64, 2}[\n          [ 1.0 ]\',\n          [ 0.75 0.25 ],\n          [ 0.75 0.25 ; 0.25 0.75 ]\n      ]\n                                        ) do sp, t, i\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],\n        mupliplier * fuel_cost[t] * thermal_generation\n    )\n    if i == 1  # wet climate state\n        setnoiseprobability!(sp, [1/6, 1/3, 0.5])\n    else       # dry climate state\n        setnoiseprobability!(sp, [0.5, 1/3, 1/6])\n    end\nend"
},

{
    "location": "tutorial/05_risk.html#Formulating-the-problem-1",
    "page": "Tutorial Five: risk",
    "title": "Formulating the problem",
    "category": "section",
    "text": "For this problem, we are going to use a convex combination of the expectation ($\\mathbb{E}$) and average value-at-risk measures (AV@R${}_{1-\\beta}$). In particular, we use AV@R at the β=0.1 quantile (i.e. the worst 10% of outcomes). This can be constructed as:risk_measure = 0.5 * Expectation() + 0.5 * AVaR(0.1)Since this is a commonly used risk measure, a slightly more computationally efficient form is EAVaR:risk_measure = EAVaR(lambda=0.5, beta=0.1)This is short-hand for lambda * Expectation() + (1-lambda) * AVaR(beta). As lambda and beta tend toward 1.0, the measure becomes more risk-neutral (i.e. less risk averse).Risk measures are set in the model using the risk_measure keyword in the SDDPModel constructor. For example, our model is now:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0,\n      markov_transition = Array{Float64, 2}[\n          [ 1.0 ]\',\n          [ 0.75 0.25 ],\n          [ 0.75 0.25 ; 0.25 0.75 ]\n      ],\n           risk_measure = EAVaR(lambda=0.5, beta=0.1)\n                                        ) do sp, t, i\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],\n        mupliplier * fuel_cost[t] * thermal_generation\n    )\n    if i == 1  # wet climate state\n        setnoiseprobability!(sp, [1/6, 1/3, 0.5])\n    else       # dry climate state\n        setnoiseprobability!(sp, [0.5, 1/3, 1/6])\n    end\nend"
},

{
    "location": "tutorial/05_risk.html#Solving-the-problem-1",
    "page": "Tutorial Five: risk",
    "title": "Solving the problem",
    "category": "section",
    "text": "Deciding when to terminate a risk-averse SDDP model is an unresolved problem in the literature. In addition to the termination methods discussed in previous tutorials, the user can also terminate the solve using the keys [CRTL]+[C].In addition, to demonstrate that we cannot use Monte Carlo simulation to estimate the upper bound of a risk-averse model, we perform a Monte Carlo simulation of the policy every two iterations. If left unchecked, this solve will not terminate as we have set terminate=false.status = solve(m,\n    simulation = MonteCarloSimulation(\n        frequency  = 2,\n        confidence = 0.95,\n        terminate  = false,\n        min        = 50,\n        step       = 50,\n        max        = 100,\n    )\n)After 7 iterations, we interrupt the solve using the [CRTL]+[C] keys. (Since this is a trivial model to solve, we had to be very quick to terminate!) The return value status is :interrupted, and the log is:-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Serial solver\n    Model:\n        Stages:         3\n        States:         1\n        Subproblems:    5\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n       21.000K         8.520K        |     1    0.0      0    0.0    0.0\n   6.857K    9.363K   12.465K -45.0  |     2    0.0     50    0.0    0.0\n        7.000K        12.465K        |     3    0.0     50    0.0    0.1\n   7.911K   11.189K   12.465K -36.5  |     4    0.0    100    0.1    0.1\n        8.000K        12.477K        |     5    0.0    100    0.1    0.1\n   7.490K   10.230K   12.477K -40.0  |     6    0.0    150    0.1    0.1\n        2.000K        12.477K        |     7    0.0    150    0.1    0.1\nWARNING: Terminating solve due to user interaction\n-------------------------------------------------------------------------------\n    Other Statistics:\n        Iterations:         7\n        Termination Status: interrupted\n==============================================================================="
},

{
    "location": "tutorial/05_risk.html#Extra-for-experts:-new-risk-measures-1",
    "page": "Tutorial Five: risk",
    "title": "Extra for experts: new risk measures",
    "category": "section",
    "text": "One of the cool features of SDDP.jl is how easy it is to create new risk measures. To illustrate this, we consider implementing the worst-case risk measure. First, we need to create a new concrete subtype of the abstract type AbstractRiskMeasure defined by SDDP.jl:\"\"\"\n    TheWorstCase()\n\nCreate an instance of the worst-case risk measures. This places all of the\nweight on the maximum outcome if minimizing, and the minimum outcome if\nmaximizing.\n\"\"\"\nstruct TheWorstCase <: SDDP.AbstractRiskMeasure endThen, we need to overload the SDDP.modifyprobability! function provided by SDDP.jl. This function takes six arguments:an instance of the risk measure (e.g. TheWorstCase());\na vector of the risk-adjusted probability distribution that the function modifies in-place;\na vector of the original probability distribution;\na vector of the observations;\nthe SDDPModel m; and\nthe JuMP subproblem sp.For example, the worst-case risk measure places all of the probability on the worst outcome:function SDDP.modifyprobability!(::TheWorstCase,\n    risk_adjusted_distribution,\n    original_distribution::Vector{Float64},\n    observations::Vector{Float64},\n    m::SDDPModel,\n    sp::JuMP.Model\n    )\n    if getsense(sp) == :Min\n        worst_index = indmax(observations)\n    else\n        worst_index = indmin(observations)\n    end\n    risk_adjusted_distribution .= 0.0\n    risk_adjusted_distribution[worst_index] = 1.0\n    return nothing\nendnote: Note\nThis implementation isn\'t a proper implementation as it assumes that the worst-case outcome has a positive probability of occurring. Accounting for this edge-case efficiently makes the implementation to verbose for this simple example.Now TheWorstCase() can be used like a risk measure defined by SDDP.jl. It is even possible to compose it with other risk measures, for example:risk_measure = 0.5 * Expectation() + 0.5 * TheWorstCase()This concludes our fifth tutorial for SDDP.jl. In the next tutorial, Tutorial Six: cut selection, we introduce cut selection."
},

{
    "location": "tutorial/06_cut_selection.html#",
    "page": "Tutorial Six: cut selection",
    "title": "Tutorial Six: cut selection",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/06_cut_selection.html#Tutorial-Six:-cut-selection-1",
    "page": "Tutorial Six: cut selection",
    "title": "Tutorial Six: cut selection",
    "category": "section",
    "text": "Over the previous five tutorials, we formulated a simple risk-averse hydrothermal scheduling problem. Now, in this tutorial, we improve the computational performance of the solution process using cut selection.As the SDDP algorithm progresses, many cuts (typically thousands) are added to each subproblem. This increases the computational effort required to solve each subproblem. In addition, many of the cuts created early in the solution process may be redundant once additional cuts are added. This issue has spurred a vein of research into heuristics for choosing cuts to keep or remove (this is referred to as cut selection). To facilitate the development of cut selection heuristics, SDDP.jl features the concept of a cut oracle associated with each subproblem. A cut oracle has two jobs: it should store the complete list of cuts created for that subproblem; and when asked, it should provide a subset of those cuts to retain in the subproblem.So far, only the Level One cut selection method is implemented. It can be constructed using LevelOneCutOracle. Cut selection can be added to a model using the cut_oracle keyword argument to the SDDPModel constructor.Therefore, our risk-averse multistage stochastic optimization problem with cut selection can be formulated in SDDP.jl as:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0,\n      markov_transition = Array{Float64, 2}[\n          [ 1.0 ]\',\n          [ 0.75 0.25 ],\n          [ 0.75 0.25 ; 0.25 0.75 ]\n      ],\n           risk_measure = EAVaR(lambda=0.5, beta=0.1),\n             cut_oracle = LevelOneCutOracle()\n                                        ) do sp, t, i\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],\n        mupliplier * fuel_cost[t] * thermal_generation\n    )\n    if i == 1  # wet climate state\n        setnoiseprobability!(sp, [1/6, 1/3, 0.5])\n    else       # dry climate state\n        setnoiseprobability!(sp, [0.5, 1/3, 1/6])\n    end\nend"
},

{
    "location": "tutorial/06_cut_selection.html#Solving-the-problem-1",
    "page": "Tutorial Six: cut selection",
    "title": "Solving the problem",
    "category": "section",
    "text": "info: Info\nThis will change once JuMP 0.19 lands.Due to the design of JuMP, we are unable to delete cuts from the model. Therefore, selecting a subset of cuts involves rebuilding the subproblems from scratch. The user can control the frequency by which the cuts are selected and the subproblems rebuilt with the cut_selection_frequency keyword argument of the solve method. Frequent cut selection (i.e. when cut_selection_frequency is small) reduces the size of the subproblems that are solved but incurs the overhead of rebuilding the subproblems. However, infrequent cut selection (i.e. when cut_selection_frequency is large) allows the subproblems to grow large (by adding many constraints), leading to an increase in the solve time of individual subproblems. Ultimately, this will be a model-specific trade-off. As a rule of thumb, simpler (i.e. few variables and constraints) models benefit from more frequent cut selection compared with complicated (i.e. many variables and constraints) models.status = solve(m,\n    iteration_limit         = 10,\n    cut_selection_frequency = 5\n)We can get the subproblem from a SDDPModel using SDDP.getsubproblem. For example, the subproblem in the first stage and first Markov state is:sp = SDDP.getsubproblem(m, 1, 1)Then, given a subproblem sp, we can get the cut oracle as follows:oracle = SDDP.cutoracle(sp)Finally, we can query the list of valid cuts in the oracle using SDDP.validcuts:julia> SDDP.validcuts(oracle)\n1-element Array{Cut, 1}:\n SDDP.Cut(29964.8, [-108.333])Whereas we can query all of the cuts in the oracle using SDDP.allcuts:julia> SDDP.allcuts(oracle)\n10-element Array{Cut, 1}:\n SDDP.Cut(29548.2, [-108.229])\n SDDP.Cut(29548.2, [-108.229])\n SDDP.Cut(29964.8, [-108.333])\n SDDP.Cut(29964.8, [-108.333])\n SDDP.Cut(29964.8, [-108.333])\n ⋮So, despite performing 10 SDDP iterations, only one cut is needed to approximate the cost-to-go function! A similar result can be found in the wet Markov state in the second stage:julia> SDDP.validcuts(SDDP.cutoracle(SDDP.getsubproblem(m, 2 1)))\n3-element Array{Cut, 1}:\n SDDP.Cut(20625.0, [-162.5])\n SDDP.Cut(16875.0, [-112.5])\n SDDP.Cut(19375.0, [-137.5])"
},

{
    "location": "tutorial/06_cut_selection.html#Extra-for-experts:-new-cut-selection-heurisitics-1",
    "page": "Tutorial Six: cut selection",
    "title": "Extra for experts: new cut selection heurisitics",
    "category": "section",
    "text": "One of the cool features of SDDP.jl is how easy it is to create new cut selection heuristics.To illustrate this, we implement the Last Cuts strategy. This heuristic selects the last N discovered cuts to keep in each subproblem. First, we need to create a new concrete subtype of the abstract type AbstractCutOracle defined by SDDP.jl:\"\"\"\n    LastCuts(N::Int)\n\nCreate a cut oracle that keeps the last `N` discovered cuts.\n\"\"\"\nstruct LastCuts <: SDDP.AbstractCutOracle\n    cuts::Vector{SDDP.Cut}\n    N::Int\n    LastCuts(N::Int) = new(SDDP.Cut[], N)\nendThen, we need to overload three methods: SDDP.storecut!, SDDP.validcuts, and SDDP.allcuts.SDDP.storecut takes four arguments. The first is an instance of the cut oracle. The second and third arguments are the SDDPModel m and the JuMP subproblem sp. The fourth argument is the cut itself. For our example, we store all the cuts that have been discovered:function SDDP.storecut!(oracle::LastCuts, m::SDDPModel, sp::JuMP.Model, cut::SDDP.Cut)\n    push!(oracle.cuts, cut)\nendSDDP.validcuts returns a list of all of the cuts that are valid cuts to keep in the subproblem. The LastCuts oracle returns the most recent N discovered cuts:function SDDP.validcuts(oracle::LastCuts)\n    oracle.cuts[max(1, end - oracle.N + 1):end]\nendFinally, SDDP.allcuts returns a list of all of the cuts that have been discovered:function SDDP.allcuts(oracle::LastCuts)\n    oracle.cuts\nendNow, the cut oracle LastCuts(500) can be used just like any other cut oracle defined by SDDP.jl.This concludes our sixth tutorial for SDDP.jl. In our next tutorial, Tutorial Seven: plotting, we discuss some of the plotting utilities of SDDP.jl"
},

{
    "location": "tutorial/07_plotting.html#",
    "page": "Tutorial Seven: plotting",
    "title": "Tutorial Seven: plotting",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/07_plotting.html#Tutorial-Seven:-plotting-1",
    "page": "Tutorial Seven: plotting",
    "title": "Tutorial Seven: plotting",
    "category": "section",
    "text": "In our previous tutorial, we have formulated, solved, and simulated a multistage stochastic optimization problem. However, we haven\'t really investigated what the solution looks like. Luckily, SDDP.jl includes a number of plotting tools to help us do that. In this tutorial, we explain the tools and make some pretty pictures.First, recall from Tutorial Six: cut selection that our model is:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0,\n      markov_transition = Array{Float64, 2}[\n          [ 1.0 ]\',\n          [ 0.75 0.25 ],\n          [ 0.75 0.25 ; 0.25 0.75 ]\n      ],\n           risk_measure = EAVaR(lambda=0.5, beta=0.1),\n             cut_oracle = LevelOneCutOracle()\n                                        ) do sp, t, i\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],\n        mupliplier * fuel_cost[t] * thermal_generation\n    )\n    if i == 1  # wet climate state\n        setnoiseprobability!(sp, [1/6, 1/3, 0.5])\n    else       # dry climate state\n        setnoiseprobability!(sp, [0.5, 1/3, 1/6])\n    end\nendWe\'re going to solve this for 20 iterations and then simulate 100 Monte Carlo realizations of the solution.status = solve(m, iteration_limit = 20)\nsimulation_result = simulate(m,\n    100,\n    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill]\n)"
},

{
    "location": "tutorial/07_plotting.html#Plotting-the-simulated-trajectories-1",
    "page": "Tutorial Seven: plotting",
    "title": "Plotting the simulated trajectories",
    "category": "section",
    "text": "inflows = [0.0, 50.0, 100.0]\nplt = SDDP.newplot()\nSDDP.addplot!(plt,\n    1:100, 1:3,\n    (i, t)->simulation_result[i][:thermal_generation][t],\n    title  = \"Thermal Generation\",\n    ylabel = \"MWh\"\n)\nSDDP.addplot!(plt,\n    1:100, 1:3,\n    (i, t)->inflows[simulation_result[i][:noise][t]],\n    title  = \"Inflows\",\n    ylabel = \"MWh\"\n)\nSDDP.show(plt)This should open a plot window with a plot that looks like:(Image: single_trajectory)Using the mouse, you can highlight individual trajectories by hovering over them. This makes it possible to visualize a single trajectory across multiple dimensions. In the above plot, we are hovering over the highest thermal generation trajectory. As could be expected, this occurs when the inflow is 0 in every stage.If you click on the plot, then trajectories that are close to the mouse pointer are shown darker and those further away are shown lighter. In the following image, we clicked on the high thermal generation point in the first stage. This shows that thermal generation is high when inflows are low.(Image: clicked_trajectory)"
},

{
    "location": "tutorial/07_plotting.html#Publication-quality-plots-1",
    "page": "Tutorial Seven: plotting",
    "title": "Publication quality plots",
    "category": "section",
    "text": "Instead of the interactive Javascript plots, you can also create some publication ready plots using the SDDP.publicationplot function.note: Note\nYou need to install the Plots.jl package for this to work. We used the GR backend (gr()), but any Plots.jl backend should work.This function implements a plot recipe to create ribbon plots of each variable against the stages. The first argument is the vector of simulation dictionaries and the second argument is the dictionary key that you want to plot. Standard Plots.jl keyword arguments such as title and xlabel can be used to modify the look of each plot. By default, the plot displays ribbons of the 0-100, 10-90, and 25-75 percentiles. The dark, solid line in the middle is the median (i.e. 50\'th percentile).using Plots\ngr()\nplot(\n    SDDP.publicationplot(simulation_result, :outgoing_volume, title=\"Volume\"),\n    SDDP.publicationplot(simulation_result, :thermal_generation, title=\"Thermal Generation\"),\n    SDDP.publicationplot(simulation_result, :hydro_generation, title=\"Hydro Generation\"),\n    SDDP.publicationplot(simulation_result, :hydro_spill, title=\"Hydro Spill\"),\n    layout        = (2,2),\n    size          = (1000, 600),\n    titlefont     = Plots.font(\"times\", 14),\n    guidefont     = Plots.font(\"times\", 14),\n    tickfont      = Plots.font(\"times\", 14),\n    bottom_margin = 7.5Plots.mm,\n    left_margin   = 5Plots.mm,\n    xlabel        = \"Stage\\n\",\n    xticks        = [1,2,3]\n)This should open a plot window with a plot that looks like:(Image: publication plot)You can save this to a PDF using Plots.jlsavefig(\"my_picture.pdf\")"
},

{
    "location": "tutorial/07_plotting.html#Plotting-the-value-function-1",
    "page": "Tutorial Seven: plotting",
    "title": "Plotting the value function",
    "category": "section",
    "text": "It is often of interest to examine the shape of the cost-to-go function. SDDP.jl facilitates this with SDDP.plotvaluefunction. It requires at least four inputs. The first is the SDDPModel m, the second is the index of the stage, and the third is the index of the Markov state. The remaining arguments define a set of discretized points at which to evaluate the cost-to-go function. The fourth argument gives the set of points at which to evaluate the first state variable, the fifth argument gives the points for the second state variable (if one exists), and so on.For our example, there is only one state variable, so to visualize it over the range of 0:200 in the wet Markov state in the second stage, we call:SDDP.plotvaluefunction(m, 2, 1, 0.0:200.0; label1=\"Volume\")This should open a browser window with a plot that looks like: (Image: value function)SDDP.plotvaluefunction can also be used to visualize problems with multiple state dimensions. For example, if we had three reservoirs, we can fix the second reservoir to 100.0 units and then visualize the cost-to-go surface with respect to the two free dimensions:SDDP.plotvaluefunction(m, 2, 1,\n    0.0:1:200.0,\n    100.0,\n    0.0:1:200.0;\n    label1=\"First State\",\n    label2=\"Third State\"\n)That concludes our seventh tutorial for SDDP.jl. In our next tutorial, Tutorial Eight: odds and ends, we  discuss some odds and ends relating to SDDP now that we have a basic  understanding of it\'s functionality."
},

{
    "location": "tutorial/08_odds_and_ends.html#",
    "page": "Tutorial Eight: odds and ends",
    "title": "Tutorial Eight: odds and ends",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/08_odds_and_ends.html#Tutorial-Eight:-odds-and-ends-1",
    "page": "Tutorial Eight: odds and ends",
    "title": "Tutorial Eight: odds and ends",
    "category": "section",
    "text": "In our previous tutorials, we have discussed the formulation, solution, and visualization of a multistage stochastic optimization problem. By now you should have all the tools necessary to solve your own problems using SDDP.jl. However, there are a few odds and ends that we need to address."
},

{
    "location": "tutorial/08_odds_and_ends.html#The-@state-macro-1",
    "page": "Tutorial Eight: odds and ends",
    "title": "The @state macro",
    "category": "section",
    "text": "In Tutorial One: first steps, we introduced a single state variable. However, much more complex syntax is supported. This syntax hijacks JuMP\'s @variable macro syntax, so it will be familiar to users of JuMP, but may look unusual at first glance.RESERVOIRS = [ :upper, :middle, :lower ]\nV0 = Dict(\n    :upper  = 10.0,\n    :middle =  5.0,\n    :lower  =  2.5\n)\n@state(m, volume[reservoir=RESERVOIRS], initial_volume == V0[reservoir])This will create JuMP variables that can be accessed as volume[:upper] and initial_volume[:lower]. There is also @states, an analogue to JuMP\'s @variables macro. It can be used as follows:@states(sp, begin\n    x >= 0,   x0==1\n    y[i=1:2], y0==i\nend)"
},

{
    "location": "tutorial/08_odds_and_ends.html#Multiple-RHS-noise-constraints-1",
    "page": "Tutorial Eight: odds and ends",
    "title": "Multiple RHS noise constraints",
    "category": "section",
    "text": "In Tutorial Two: RHS noise, we added a single constraint with noise in the RHS term; however, you probably want to add many of these constraints. There are two ways to do this. First, you can just add multiple calls like:@rhsnoise(sp, w=[1,2,3], x <= w)\n@rhsnoise(sp, w=[4,5,6], y >= w)Multiple calls to @rhsnoise, they are _not_ sampled independently. Instead, in the example above, there will be three possible realizations: (1,4), (2,5), and (3,6). Therefore, if you have multiple calls to @rhsnoise, they must have the same number of elements in their sample space. For example, the following will not work:@rhsnoise(sp, w=[1,2,3], x <= w)    # 3 elements\n@rhsnoise(sp, w=[4,5,6,7], y >= w)  # 4 elementsAnother option is to use the @rhsnoises macro. It is very similar to @states and JuMP\'s @consraints macro:@rhsnoises(sp, w=[1,2,3], begin\n    x <= w\n    y >= w + 3\nend)"
},

{
    "location": "tutorial/08_odds_and_ends.html#if-statements-in-more-detail-1",
    "page": "Tutorial Eight: odds and ends",
    "title": "if statements in more detail",
    "category": "section",
    "text": "In Tutorial Four: Markovian policy graphs, we used an if statement to control the call to setnoiseprobability! depending on the Markov state. This can be used much more generally. For example:m = SDDPModel(\n    #...arguments omitted ...\n        ) do sp, t, i\n    @state(m, x>=0, x0==0)\n    if t == 1\n        @stageobjective(m, x)\n    else\n        if i == 1\n            @variable(m, u >= 1)\n            @constraint(m, x0 + u == x)\n        else\n            @rhsnoise(m, w=[1,2], x0 + w == x)\n        end\n        @stageobjective(m, x0 + x)\n    end\nendYou could, of course, do the following as well:function first_stage(sp::JuMP.Model, x0, x)\n    @stageobjective(m, x)\nend\nfunction second_stage(sp::JuMP.Model, x0, x)\n    if i == 1\n        @variable(m, u >= 1)\n        @constraint(m, x0 + u == x)\n    else\n        @rhsnoise(m, w=[1,2], x0 + w == x)\n    end\n    @stageobjective(m, x0 + x)\nend\nm = SDDPModel(\n    #...arguments omitted ...\n        ) do sp, t, i\n    @state(m, x>=0, x0==0)\n    if t == 1\n        first_stage(m, x0, x)\n    else\n        second_stage(m, x0, x)\n    end\nend"
},

{
    "location": "tutorial/08_odds_and_ends.html#Stage-dependent-risk-measures-1",
    "page": "Tutorial Eight: odds and ends",
    "title": "Stage-dependent risk measures",
    "category": "section",
    "text": "In Tutorial Five: risk, we discussed adding risk measures to our model. We used the same risk measure for every stage. However, often you may want to have time-varying risk measures. This can be accomplished very easily: instead of passing a risk measure to the risk_measure keyword argument, we pass a function that takes two inputs: t::Int, the index of the stage; and i::Int, the index of the Markov state. For example:function build_my_risk_measure(t::Int, i::Int)\n    if t == 1\n        return Expectation()\n    elseif i == 1\n        return AVaR(0.5)\n    else\n        return 0.5 * Expectation() + 0.5 * WorstCase()\n    end\nend\n\nm = SDDPModel(\n    # ... arguments omitted ...\n    risk_measure = build_my_risk_measure\n                                ) do sp, t, i\n    # ... formulation omitted ...\nend"
},

{
    "location": "tutorial/08_odds_and_ends.html#Reading-and-writing-cuts-to-file-1",
    "page": "Tutorial Eight: odds and ends",
    "title": "Reading and writing cuts to file",
    "category": "section",
    "text": "It\'s possible to save the cuts that are discovered during a solve to a file so that they can be read back in or used for other analysis. This can be done using the cut_output_file to solve:m = build_model()\nSDDP.solve(m,\n    iteration_limit = 10,\n    cut_output_file = \"cuts.csv\"\n)cuts.csv is a csv file with N+3 columns, where N is the number of state variables in the model. The columns are: (1) the index of the stage; (2) the index of the Markov state; (3) the intercept; and (4+), the coefficients of the cuts for each of the state variables.This cut file can be read back into a model using loadcuts!:m2 = build_model()\nloadcuts!(m2, \"cuts.csv\")Another option is to use the SDDP.writecuts! method after the model has been solved:m = build_model()\nSDDP.solve(m, iteration_limit = 10)\nSDDP.writecuts!(\"cuts.csv\", m)"
},

{
    "location": "tutorial/08_odds_and_ends.html#Timer-outputs-1",
    "page": "Tutorial Eight: odds and ends",
    "title": "Timer outputs",
    "category": "section",
    "text": "SDDP.jl embeds the great TimerOutputs.jl package to help profile where time is spent during the solution process. You can print the timing statistics by setting print_level=2 in a solve call. This produces a log like:-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Serial solver\n    Model:\n        Stages:         3\n        States:         1\n        Subproblems:    5\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n       27.000K         5.818K        |     1    0.0      0    0.0    0.0\n        2.000K         6.952K        |     2    0.0      0    0.0    0.0\n        2.000K         6.952K        |     3    0.0      0    0.0    0.0\n       11.000K         7.135K        |     4    0.0      0    0.0    0.0\n        2.000K         7.135K        |     5    0.0      0    0.0    0.0\n        2.000K         7.135K        |     6    0.0      0    0.0    0.0\n        5.000K         7.135K        |     7    0.0      0    0.0    0.0\n        2.000K         7.135K        |     8    0.0      0    0.0    0.0\n        5.000K         7.135K        |     9    0.0      0    0.0    0.0\n       12.500K         7.135K        |    10    0.0      0    0.0    0.0\n-------------------------------------------------------------------------------\n ─────────────────────────────────────────────────────────────────────────────────\n         Timing statistics                Time                   Allocations\n                                  ──────────────────────   ───────────────────────\n         Tot / % measured:            40.8ms / 98.0%           0.97MiB / 100%\n\n Section                  ncalls     time   %tot     avg     alloc   %tot      avg\n ─────────────────────────────────────────────────────────────────────────────────\n Solve                         1   40.0ms   100%  40.0ms   0.96MiB  100%   0.96MiB\n   Iteration Phase            10   38.3ms  95.6%  3.83ms    950KiB  96.3%  95.0KiB\n     Backward Pass            10   30.5ms  76.3%  3.05ms    784KiB  79.4%  78.4KiB\n       JuMP.solve            150   23.8ms  59.4%   159μs    423KiB  42.9%  2.82KiB\n         optimize!           150   19.7ms  49.2%   131μs         -  0.00%        -\n         prep JuMP model     150   2.69ms  6.72%  17.9μs    193KiB  19.6%  1.29KiB\n         getsolution         150   1.09ms  2.71%  7.23μs    209KiB  21.1%  1.39KiB\n       Cut addition           30   1.12ms  2.81%  37.4μs   69.1KiB  7.01%  2.30KiB\n         risk measure         30   27.5μs  0.07%   917ns   8.91KiB  0.90%        -\n     Forward Pass             10   7.68ms  19.2%   768μs    163KiB  16.5%  16.3KiB\n       JuMP.solve             30   5.89ms  14.7%   196μs   93.4KiB  9.47%  3.11KiB\n         optimize!            30   4.59ms  11.5%   153μs         -  0.00%        -\n         prep JuMP model      30    909μs  2.27%  30.3μs   45.6KiB  4.62%  1.52KiB\n         getsolution          30    303μs  0.76%  10.1μs   41.6KiB  4.21%  1.39KiB\n ─────────────────────────────────────────────────────────────────────────────────\n    Other Statistics:\n        Iterations:         10\n        Termination Status: iteration_limit\n==============================================================================="
},

{
    "location": "tutorial/08_odds_and_ends.html#Discount-factors-1",
    "page": "Tutorial Eight: odds and ends",
    "title": "Discount factors",
    "category": "section",
    "text": "If a model has a large number of stages, a common modelling trick is to add a discount factor to the cost-to-go function. The stage-objective can then be written as c^top + rho theta, where \\\\theta is the cost-to-go variable. However, since SDDP.jl does not expose the cost-to-go variable to the user, this modelling trick must be accomplished as follows:SDDPModel() do sp, t\n    ρ = 0.95\n    # ... lines omitted ...\n    @stageobjective(sp, ρ^(t - 1) * dot(c, x))\nendThat concludes our eighth tutorial for SDDP.jl. In our next tutorial, Tutorial Nine: nonlinear models, we discuss how SDDP.jl can be used to solve problems that have nonlinear transition functions."
},

{
    "location": "tutorial/09_nonlinear.html#",
    "page": "Tutorial Nine: nonlinear models",
    "title": "Tutorial Nine: nonlinear models",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/09_nonlinear.html#Tutorial-Nine:-nonlinear-models-1",
    "page": "Tutorial Nine: nonlinear models",
    "title": "Tutorial Nine: nonlinear models",
    "category": "section",
    "text": "In our previous tutorials, we formulated a linear version of the hydrothermal scheduling problem. To do so, we had to make a large assumption, namely, that the inflows were stagewise-independent. In Tutorial Four: Markovian policy graphs, we improved upon this slightly using a Markov chain with persistent climate states. However, another commonly used model is to assume that the inflows follow a log auto-regressive process. In this tutorial, we explain how to implement this in SDDP.jl.As a starting point, we use a model that is very similar to the hydrothermal scheduling problem we formulated in Tutorial Two: RHS noise. In that model, we assumed that there the a constraint:@rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n    outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n)We assumed that the inflow term was stagewise-independent. In this tutorial, we assume that the inflow term can be modelled like:log(inflow[t]) = log(inflow[t-1]) + log(noise),where noise is drawn from [0.9, 1.0, 1.1] with uniform probability. Since we value of the inflow in the previous stage affects the value in the current stage, it needs to be added as a state variable:@state(sp, current_inflow >= 1, previous_inflow == 50)Note that we place a small, positive lower bound on the current_inflow to prevent the solver calling the undefined log(0.0).Now we want to add the log transition constraint. We can do this using JuMP\'s @NLconstraint macro to add this constraint. However, SDDP.jl only supports noise terms in the right hand-side of a linear constraint. We can overcome this limitation by introducing a dummy variable (note that it also has a small, positive lower bound to avoid log(0.0)):@variable(sp, noise >= 0.1)\n@rhsnoise(sp, ω = [0.9, 1.0, 1.1], noise == ω)\n@NLconstraint(sp, log(current_inflow) == log(previous_inflow) + log(noise))Finally, we\'re using JuMP\'s nonlinear functionality, we need to choose an appropriate solver. We choose to use the COIN solver Ipopt. Therefore, our final model is:using JuMP, SDDP, Ipopt\nm = SDDPModel(\n    stages          = 3,\n    sense           = :Min,\n    solver          = IpoptSolver(print_level=0),\n    objective_bound = 0.0\n                ) do sp, t\n    @states(sp, begin\n        200 >= outgoing_volume >= 0, incoming_volume == 200\n               current_inflow  >= 1, previous_inflow ==  50\n    end)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n        inflow_noise_term  >= 0.1\n    end)\n    @constraints(sp, begin\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == current_inflow\n        hydro_generation + thermal_generation >= 150.0\n    end)\n\n    @rhsnoise(sp, ω = [0.9, 1.0, 1.1], inflow_noise_term == ω)\n    @NLconstraint(sp, log(current_inflow) == log(previous_inflow) + log(inflow_noise_term))\n\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, fuel_cost[t] * thermal_generation )\nendThis problem can be solved just like any other SDDP model:status = solve(m, iteration_limit = 10)\nsimulation_result = simulate(m, 1, [:inflow′])Then, we can check that the inflows do indeed follow a log auto-regressive process.julia> simulaton_result[1][:inflow′]\n 55.0\n 60.5\n 54.45That concludes our ninth tutorial for SDDP.jl. In our next tutorial, Tutorial Ten: parallelism, we explain how to solve SDDP models in parallel."
},

{
    "location": "tutorial/10_parallel.html#",
    "page": "Tutorial Ten: parallelism",
    "title": "Tutorial Ten: parallelism",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/10_parallel.html#Tutorial-Ten:-parallelism-1",
    "page": "Tutorial Ten: parallelism",
    "title": "Tutorial Ten: parallelism",
    "category": "section",
    "text": "The SDDP algorithm is highly parallelizable. In SDDP.jl, we chose to implement an approach that minimizes inter-process communication. We call this asynchronous SDDP.In our implementation, one process is designated the master process and the remainder are designated as slaves. Each slave receives a full copy of the SDDP model and is set to work, performing iterations. At the end of each iteration, the slave passes the master the cuts it discovered during the iteration and receives any new cuts discovered by other slaves. The slave also queries the master as to whether it should terminate, perform another iteration, or perform a simulation. If the master requests a simulation (for example, to calculate a confidence interval in order to test for convergence), the slave returns the objective value of the simulation rather than a new set of cuts.In this tutorial, we explain how to use the asynchronous solve feature of SDDP.jl.First, we need to add some extra processes to Julia. This can be done two ways. We can start Julia using julia -p N, where N is the number of worker processes to add, or we can use the addprocs(N) function while Julia is running. The main process that co-ordinates everything is called the master process, and the remote  processes are called workers. For example, to add two workers, we run:julia> addprocs(2)\n2-element Array{Int64, 1}:\n 2\n 3One way of running a command on all of the processes is to use the @everywhere macro. We can also use the myid() function to query the index of the process. For example:julia> @everywhere println(\"Called from: $(myid())\")\nCalled from: 1\n        From worker 3:  Called from: 3\n        From worker 2:  Called from: 2Note that the order we receive things from the worker processes is not deterministic.Now what we need to do is to initialize the random number generator on each process. However, we need to be careful not to use the same seed on each process or each process will perform identical SDDP iterations!julia> @everywhere srand(123 * myid())This will set srand(123) on process 1, srand(246) on process 2, and so on.!!!note     If you are using any data or function in the subproblem definition, these     need to be copied to every process. The easiest way to do this is to place     everything in a file and then run @everywhere include(\"path/to/my/file\").Recall from Tutorial Four: Markovian policy graphs that our model is:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0,\n      markov_transition = Array{Float64, 2}[\n          [ 1.0 ]\',\n          [ 0.75 0.25 ],\n          [ 0.75 0.25 ; 0.25 0.75 ]\n      ]\n                                        ) do sp, t, i\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],\n        mupliplier * fuel_cost[t] * thermal_generation\n    )\n    if i == 1  # wet climate state\n        setnoiseprobability!(sp, [1/6, 1/3, 0.5])\n    else       # dry climate state\n        setnoiseprobability!(sp, [0.5, 1/3, 1/6])\n    end\nendWe can solve this using SDDP.jl\'s asynchronous feature by passing an instance of Asynchronous to the solve_type keyword in solve:status = solve(m,\n    iteration_limit = 10,\n    solve_type      = Asynchronous()\n)If you have multiple processes, SDDP.jl will detect this and choose asynchronous by default. You can force the serial solution by passing an instance of Serial to solve_type.The log is:-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Asynchronous solver with 2 slave processors\n    Model:\n        Stages:         3\n        States:         1\n        Subproblems:    5\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n       27.000K         5.818K        |     2    0.0      0    0.0    0.0\n        0.000          6.198K        |     1    0.0      0    0.0    0.0\n        9.000K         6.952K        |     3    0.0      0    0.0    0.0\n        2.000K         7.135K        |     4    0.0      0    0.0    0.0\n        2.000K         7.135K        |     5    0.0      0    0.0    0.0\n        5.000K         7.135K        |     6    0.0      0    0.0    0.0\n        5.000K         7.135K        |     7    0.0      0    0.0    0.0\n        2.000K         7.135K        |     8    0.0      0    0.0    0.1\n       24.000K         7.135K        |     9    0.0      0    0.0    0.1\n       20.000K         7.135K        |    10    0.1      0    0.0    0.1\n-------------------------------------------------------------------------------\n    Other Statistics:\n        Iterations:         10\n        Termination Status: iteration_limit\n===============================================================================Note that the order of the Cut # column is not sequential because they are numbered in order of when they were created.This concludes our tenth tutorial on SDDP.jl. In the next tutorial, Tutorial Eleven: distributionally robust SDDP, we will learn how distributionally robust optimization can be incorporated in SDDP.jl."
},

{
    "location": "tutorial/11_DRO.html#",
    "page": "Tutorial Eleven: distributionally robust SDDP",
    "title": "Tutorial Eleven: distributionally robust SDDP",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/11_DRO.html#Tutorial-Eleven:-distributionally-robust-SDDP-1",
    "page": "Tutorial Eleven: distributionally robust SDDP",
    "title": "Tutorial Eleven: distributionally robust SDDP",
    "category": "section",
    "text": "In Tutorial Five: risk, we saw how risk measures can be used within our model. In this tutorial we will learn how to incorporate a distributionally robust optimization approach to SDDP in SDDP.jl.Distributionally robust optimization (DRO) is a modelling approach for optimization under uncertainty. In our setup, DRO is equivalent to a coherent risk measure, and we can apply DRO by using the risk_measure keyword we saw previously."
},

{
    "location": "tutorial/11_DRO.html#A-little-motivation-on-the-concept-1",
    "page": "Tutorial Eleven: distributionally robust SDDP",
    "title": "A little motivation on the concept",
    "category": "section",
    "text": "When we build a policy using SDDP, we use a model to represent uncertain parameters. When we come to use or evaluate our policy, the realized scenarios may not actually behave in the way we modelled our uncertainty. For example, the hydrothermal scheduling model from the previous tutorials assumed that inflows are independent between stages. However, a real sequence of inflows is likely to exhibit correlation between stages. Furthermore, we may wish to hold out a set of historical inflow sequences other than those included while generating the policy, and evaluate the performance of the policy based on its performance in the held out set. The held out set may also correspond to inflows that are not stagewise independent. With a distributionally robust approach, we avoid assuming an explicit model on the probabilities of the scenarios we consider. Instead, each time we come to add a cut, we assume that the probabilities associated with each noise are the worst case probabilities possible (with respect to our objective) within some ambiguity set.The implementation of distributionally robust SDDP here comes from the paper: A.B. Philpott, V.L. de Matos, L. Kapelevich (2018): Distributionally Robust SDDP, Computational Management Science, (link) where the details of the approach are described."
},

{
    "location": "tutorial/11_DRO.html#Formulating-the-problem-1",
    "page": "Tutorial Eleven: distributionally robust SDDP",
    "title": "Formulating the problem",
    "category": "section",
    "text": "In Tutorial Two: RHS noise, we formulated a hydrothermal scheduling problem with uncertainty in the right hand-side of the constraints. In the following, we present a similar model; however, we have set the inflow in the first stage to deterministic and equal to 50.0. The model is:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0\n                                        ) do sp, t\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    if t == 1\n        @constraint(sp, outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == 50.0)\n    else\n        @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n            outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n        )\n      setnoiseprobability!(sp, [1/3, 1/3, 1/3])\n    end\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, fuel_cost[t] * thermal_generation )\nendThis model assumed that the probability of each inflow is equally likely, with probability 1/3. The lower bound we converged to was 8.33k.To describe a distributionally robust version of this problem, we need to choose a radius of uncertainty. This is the maximum distance around the default probability vector ([1/3, 1/3, 1/3]) we will consider in our ambiguity set. For a problem with S noises, this radius should be less than $\\sqrt{(S-1)/S}$ (which would be the same as TheWorstCase() from Tutorial Five: risk).note: Note\nWe currently assume our uncertainty set is a ball centered around the probability vector assigning equal probabilities to all noises. We could generalize the algorithm to alter the center of the ball, but this is not currently implemented (see this issue).Suppose, for example, we choose the radius of uncertainty to be 1/6. Intuitively, this means we can decrease the risk-adjusted probability of each noise by almost 1/6 at the most. We can implement this by inserting a DRO(1/6) object for the keyword argument risk_measure we saw earlier.This gives the new model:m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0,\n           risk_measure = DRO(1/6)\n                                        ) do sp, t\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    if t == 1\n        @constraint(sp, outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == 50.0)\n    else\n        @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n            outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n        )\n      setnoiseprobability!(sp, [1/3, 1/3, 1/3])\n    end\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, fuel_cost[t] * thermal_generation )\nend"
},

{
    "location": "tutorial/11_DRO.html#Solving-the-problem-1",
    "page": "Tutorial Eleven: distributionally robust SDDP",
    "title": "Solving the problem",
    "category": "section",
    "text": "We can solve the above problem, terminating at our choice of iteration_limit. For example,solve(m, iteration_limit = 10)gives the following output log:-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Serial solver\n    Model:\n        Stages:         3\n        States:         1\n        Subproblems:    3\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n       10.000K        10.023K        |     1    0.0      0    0.0    0.0\n        5.000K        10.023K        |     2    0.0      0    0.0    0.0\n       12.500K        10.023K        |     3    0.0      0    0.0    0.0\n       10.000K        10.023K        |     4    0.0      0    0.0    0.0\n        5.000K        10.023K        |     5    0.0      0    0.0    0.0\n        5.000K        10.023K        |     6    0.0      0    0.0    0.0\n       17.500K        10.023K        |     7    0.0      0    0.0    0.0\n       10.000K        10.023K        |     8    0.0      0    0.0    0.0\n       12.500K        10.023K        |     9    0.0      0    0.0    0.0\n        5.000K        10.023K        |    10    0.0      0    0.0    0.0\n-------------------------------------------------------------------------------\n    Other Statistics:\n        Iterations:         10\n        Termination Status: iteration_limit\n===============================================================================We have converged to a lower bound of roughly \\$10.023k. One can check that this is a little lower than the bound from the worst case measure, which is \\$15k, but greater than the lower bound using the expectation risk measure of \\$8.33k.To run a sanity check, let us set the radius to be sufficiently large in order to match WorstCase(). Since we have S=3 scenarios, any radius larger than $\\sqrt{(3-1)/3} = \\sqrt(2/3)$ will do. (This is large enough to move from the default probability vector [0.33, 0.33, 0.33] to the worst case probability vector [1.0, 0.0, 0.0].)m = SDDPModel(\n                  sense = :Min,\n                 stages = 3,\n                 solver = ClpSolver(),\n        objective_bound = 0.0,\n           risk_measure = DRO(sqrt(2/3))\n                                        ) do sp, t\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    if t == 1\n        @constraint(sp, outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == 50.0)\n    else\n        @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n            outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n        )\n      setnoiseprobability!(sp, [1/3, 1/3, 1/3])\n    end\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    fuel_cost = [50.0, 100.0, 150.0]\n    @stageobjective(sp, fuel_cost[t] * thermal_generation )\nendIndeed, solving the model like we did above provides us with a lower bound of \\$15k.-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Serial solver\n    Model:\n        Stages:         3\n        States:         1\n        Subproblems:    3\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n       20.000K        12.500K        |     1    0.0      0    0.0    0.0\n        5.000K        15.000K        |     2    0.0      0    0.0    0.0\n       10.000K        15.000K        |     3    0.0      0    0.0    0.0\n       15.000K        15.000K        |     4    0.0      0    0.0    0.0\n        5.000K        15.000K        |     5    0.0      0    0.0    0.0\n        5.000K        15.000K        |     6    0.0      0    0.0    0.0\n       15.000K        15.000K        |     7    0.0      0    0.0    0.0\n       15.000K        15.000K        |     8    0.0      0    0.0    0.0\n       10.000K        15.000K        |     9    0.0      0    0.0    0.0\n       10.000K        15.000K        |    10    0.0      0    0.0    0.0\n-------------------------------------------------------------------------------\n    Other Statistics:\n        Iterations:         10\n        Termination Status: iteration_limit\n===============================================================================That concludes our eleventh tutorial for SDDP.jl. In the next tutorial, Tutorial Twelve: price interpolation, we discuss an extension of SDDP to models with stagewise-dependent objective uncertainty."
},

{
    "location": "tutorial/12_price_interpolation.html#",
    "page": "Tutorial Twelve: price interpolation",
    "title": "Tutorial Twelve: price interpolation",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/12_price_interpolation.html#Tutorial-Twelve:-price-interpolation-1",
    "page": "Tutorial Twelve: price interpolation",
    "title": "Tutorial Twelve: price interpolation",
    "category": "section",
    "text": "There are many applications in which we want to model a price process that follows some auto-regressive process. Common examples include stock prices on financial exchanges and spot-prices in energy markets. However, it is well known that these cannot be incorporated in to SDDP because they result in cost-to-go functions that are convex with respect to some state variables (e.g., the reservoir levels) and concave with respect to other state variables (e.g., the spot price in the current stage). To overcome this problem, the approach in the literature has been to discretize the price process in order to model it using a Markovian policy graph like those discussed in Tutorial Four: Markovian policy graphs.However, recent work offers a way to include stagewise-dependent objective uncertainty into the objective function of SDDP subproblems. Readers are directed to the following works for an introduction:Downward, A., Dowson, O., and Baucke, R. (2017). Stochastic dual dynamic programming with stagewise dependent objective uncertainty. Optimization Online.\nDowson, O. PhD Thesis. University of Auckland, 2018. linkIn Tutorial Five: risk, we formulated a risk-averse version of the hydrothermal scheduling problem. In this tutorial, we extend that model to the case where the fuel cost follows the log auto-regressive process: log(fuel_cost[t]) = log(fuel_cost[t-1]) + log(noise)where noise is drawn from the sample space [0.9, 1.0, 1.1] with equal probability.To model this in SDDP.jl, we can pass a DynamicPriceInterpolation object to the value_function keyword in SDDPModel. DynamicPriceInterpolation takes a number of arguments. First, we need to pass dynamics a function that takes two inputs – the value of the price state in the previous stage and an instance of the  noise – and returns value of the price state for the current stage. For example, the dynamics of our price process are:function fuel_cost_dynamics(fuel_cost, noise)\n    return noise * fuel_cost\nendWe also need to specify the distribution of the noise term. We do this by passing a DiscreteDistribution to the noise keyword. DiscreteDistribution takes two arguments: the first is a vector of realizations, and the second is a corresponding vector of probabilities. For our example, we create the noise distribution as:noise = DiscreteDistribution( [0.9, 1.0, 1.1], [1/3, 1/3, 1/3] )It is the realizations of the noise 0.9, 1.0, or 1.1 that are passed as noise to fuel_cost_dynamics.We also need to pass the value of the price state in the root node to initial_price, as well as the minimum (to min_price) and maximum (to max_price) possible values of the price state variable.Finally, we need to declare a lipschitz_constant. In each stage, the lipschitz_constant should be larger than the maximum possible absolute change in the cost-to-go function given a one-unit change  in the value of the price state variable. For example, in our model, the worst-case scenario is if we are forced to use thermal generation exclusively. In that case, we need to supply 450 MWh of energy. Therefore, a one-unit change in the value of the price-state can, at most, lead to a $450 change in the cost-to-go function. However, to be on the safe side, we choose a larger value of 1000.0.Putting all of this together, we can initialize the SDDPModel using dynamic interpolation as:m = SDDPModel(\n    # ... arguments omitted ...\n    value_function = DynamicPriceInterpolation(\n        dynamics           = fuel_cost_dynamics,\n        noise              = DiscreteDistribution([0.9, 1.0, 1.1], [1/3, 1/3, 1/3]),\n        initial_price      = 100.0,\n        min_price          =  50.0,\n        max_price          = 150.0,\n        lipschitz_constant = 1000.0\n    )\n        do sp, t\n    # ... subproblem definition ...\nendIn the subproblem definition, we use a different version of the @stageobjective function. This version takes a function that maps the price in the current stage to an expression for the stage objective. For our example, the stage-objective is: @stageobjective(sp, (fuel_cost) -> fuel_cost * thermal_generation )The next question is how to extend this notation to models in which the price process depends upon the stage or Markov state. This can be implemented in SDDP.jl following a similar approach to that we discussed in Stage-dependent risk measures. Instead of passing an instance of DynamicPriceInterpolation, we pass a function that takes two arguments – the stage t and Markov state i – and returns an instance of DynamicPriceInterpolation. For our example, if the price is deterministic in the first stage:function build_price_interpolation(t::Int, i::Int)\n    noise = if t == 1\n        DiscreteDistribution([1.0], [1.0])\n    else\n        DiscreteDistribution([0.9, 1.0, 1.1], [1/3, 1/3, 1/3])\n    end\n    DynamicPriceInterpolation(\n        dynamics           = fuel_cost_dynamics,\n        initial_price      = 100.0,\n        min_price          =  50.0,\n        max_price          = 150.0,\n        noise              = noise,\n        lipschitz_constant = 1000.0\n    )\nendPutting all this together, our model is:m = SDDPModel(\n    sense             = :Min,\n    stages            = 3,\n    solver            = ClpSolver(),\n    objective_bound   = 0.0,\n    markov_transition = Array{Float64, 2}[\n        [ 1.0 ]\',\n        [ 0.75 0.25 ],\n        [ 0.75 0.25 ; 0.25 0.75 ]\n    ],\n    risk_measure      = EAVaR(lambda=0.5, beta=0.1),\n    value_function    = build_price_interpolation\n                                        ) do sp, t, i\n    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)\n    @variables(sp, begin\n        thermal_generation >= 0\n        hydro_generation   >= 0\n        hydro_spill        >= 0\n    end)\n    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],\n        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow\n    )\n    @constraints(sp, begin\n        thermal_generation + hydro_generation == 150\n    end)\n    @stageobjective(sp, (fuel_cost) -> fuel_cost * thermal_generation )\n    if i == 1  # wet climate state\n        setnoiseprobability!(sp, [1/6, 1/3, 0.5])\n    else       # dry climate state\n        setnoiseprobability!(sp, [0.5, 1/3, 1/6])\n    end\nendNow we can solve this model as usual.status = solve(m; iteration_limit=100)When we simulate the policy, we can include the extra key :price, which records the value of the price state in each stage. For example:simulation_result = simulate(m, 100,\n    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill, :price]\n)We can check that the price follows the auto-regressive process:julia> simulation_result[1][:price]\n 100.0\n  90.0\n  99.0We can also plot the cost-to-go function using SDDP.plotvaluefunction like we discussed in Tutorial Seven: plotting:SDDP.plotvaluefunction(m, 2, 2,\n    linspace(0, 200, 50),   # the reservoir volume\n    linspace(70, 130, 50);  # the price state\n    label1=\"Volume\",\n    label2=\"Price\"\n)This will launch a browser window with the following: (Image: 3d saddle function)Note that the surface is convex with respect to the volume dimension and concave with respect to the price dimension.That concludes our twelfth tutorial for SDDP.jl. In the next tutorial, Tutorial Thirteen: constraint noise, we discuss experimental support for noise in the constraint matrix."
},

{
    "location": "tutorial/13_constraint_noise.html#",
    "page": "Tutorial Thirteen: constraint noise",
    "title": "Tutorial Thirteen: constraint noise",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/13_constraint_noise.html#Tutorial-Thirteen:-constraint-noise-1",
    "page": "Tutorial Thirteen: constraint noise",
    "title": "Tutorial Thirteen: constraint noise",
    "category": "section",
    "text": "In previous tutorials (e.g., Tutorial Three: objective noise), we explicitly noted that SDDP.jl does not support noise in the constraint matrix. However, we have recently developed a hack that works around this. It should not be considered stable, and is an imperfect solution. The new version of JuMP will enable native support for these modifications instead of the current hack.!!!note     This may break SDDP extensions such as SDDiP.jl.     It may also break some features of JuMP.If w ∈ [1,2,3] with equal probability, and we want to add the constraint:@constraint(sp, 2x + w*x <= 1)Then in the subproblem definition, we can use the un-exported SDDP.addconstraintnoise! method:wx = SDDP.addconstraintnoise!(sp, x, [1,2,3])\n@constraint(m, 2x + wx <= 1)The first line introduces a new JuMP variable named wx. Depending on the realization of the noise, this will take the value of w*x. Note that we named the variable wx, but this decision was arbitrary. SDDP.addconstraintnoise! is purposefully not exported to show that it is an experimental feature.!!!note     This requires a solver that implements the MathProgBase.changecoeffs!     method.To give a full example:using SDDP, JuMP, Gurobi\nm = SDDPModel(\n        sense           = :Min,\n        stages          = 2,\n        objective_bound = 0.0,\n        solver          = GurobiSolver(OutputFlag=0) ) do sp, t\n    @state(sp, x>=0, x0==1.0)\n    @stageobjective(sp, x)\n    # instead of @constraint(sp, x == w * x0) where w ∈ [1,2], we write:\n    r_x = SDDP.addconstraintnoise!(sp, x0, [1,2])\n    @constraint(sp, x == r_x)\nendThere are four possible scenarios that can be realized in this model:(x₁, x₂) = (1, 1) with total cost of 1+1=2\n(x₁, x₂) = (1, 2) with total cost of 1+2=3\n(x₁, x₂) = (2, 2) with total cost of 2+2=4\n(x₁, x₂) = (2, 4) with total cost of 2+4=6Therefore, the expected cost is mean([2,3,4,6]) = 3.75. If we solve the model  for five iterations, we can confirm that the model converges to this bound:julia> solve(m, iteration_limit=5)\n ... log omitted ...\njulia> getbound(m)\n 3.75The log is:-------------------------------------------------------------------------------\n                          SDDP.jl © Oscar Dowson, 2017-2018\n-------------------------------------------------------------------------------\n    Solver:\n        Serial solver\n    Model:\n        Stages:         2\n        States:         1\n        Subproblems:    2\n        Value Function: Default\n-------------------------------------------------------------------------------\n              Objective              |  Cut  Passes    Simulations   Total\n     Simulation       Bound   % Gap  |   #     Time     #    Time    Time\n-------------------------------------------------------------------------------\n        4.000          3.750         |     1    1.0      0    0.0    1.0\n        3.000          3.750         |     2    1.0      0    0.0    1.3\n        4.000          3.750         |     3    1.0      0    0.0    1.3\n        2.000          3.750         |     4    1.0      0    0.0    1.3\n        2.000          3.750         |     5    1.0      0    0.0    1.3\n-------------------------------------------------------------------------------\n    Other Statistics:\n        Iterations:         5\n        Termination Status: iteration_limit\n===============================================================================Each uncertain term in the constraint matrix must be added using the SDDP.addconstraintnoise! method. These terms can be used in conjunction with right-hand side and  objective noise. In which case there must be the same number of realizations in the constraint noise as there are in the right-hand side noise terms. For example: To give a full example:using SDDP, JuMP, Gurobi\nm = SDDPModel(\n        sense           = :Min,\n        stages          = 2,\n        objective_bound = 0.0,\n        solver          = GurobiSolver(OutputFlag=0) ) do sp, t\n    @states(sp, begin\n        x>=0, x0==1.0\n        y>=0, y0==1.0\n    end)\n    # noise in the constraint matrix\n    r_x = SDDP.addconstraintnoise!(sp, x0, [0.9,1.1])\n    r_y = SDDP.addconstraintnoise!(sp, y0, [0.8,1.2])\n    @constraint(sp, x == 0.9r_x + 0.1r_y)\n\n    # noise in the right-hand side term\n    @rhsnoise(sp, w=[0.0, 0.1], y == 0.9r_y + 0.1r_x + w)\n\n    # noise in the objective function\n    @stageobjective(sp, w=[1,2], w * x + y)\nendThis concludes our thirteenth tutorial for SDDP.jl."
},

{
    "location": "readings.html#",
    "page": "Readings",
    "title": "Readings",
    "category": "page",
    "text": ""
},

{
    "location": "readings.html#Readings-1",
    "page": "Readings",
    "title": "Readings",
    "category": "section",
    "text": "On this page, we\'ve collated a variety of papers and books we think are helpful readings that cover knowledge needed to use SDDP.jl."
},

{
    "location": "readings.html#Stochastic-Optimization-1",
    "page": "Readings",
    "title": "Stochastic Optimization",
    "category": "section",
    "text": "A general primer on Stochastic ProgrammingBirge, J.R., Louveaux, F., 2011. Introduction to Stochastic Programming,  Springer Series in Operations Research and Financial Engineering. Springer New  York, New York, NY.  doi:10.1007/978-1-4614-0237-4Some overviews of Stochastic Optimization and where it sits in relation to other fieldsPowell, W.B., 2014. Clearing the Jungle of Stochastic Optimization, in: Newman,  A.M., Leung, J., Smith, J.C., Greenberg, H.J. (Eds.), Bridging Data and  Decisions. INFORMS, pp. 109–137. link\nPowell, W.B., 2016. A Unified Framework for Optimization under Uncertainty  TutORials in Operations Research, in: Optimization Challenges in Complex,  Networked and Risky Systems. pp. 45–83. link"
},

{
    "location": "readings.html#Stochastic-Dual-Dynamic-Programming-1",
    "page": "Readings",
    "title": "Stochastic Dual Dynamic Programming",
    "category": "section",
    "text": "The original paper presenting SDDPPereira, M.V.F., Pinto, L.M.V.G., 1991. Multi-stage stochastic optimization  applied to energy planning. Mathematical Programming 52, 359–375. doi:10.1007/BF01582895The paper presenting the Markov version of SDDP implemented in this libraryPhilpott, A.B., de Matos, V.L., 2012. Dynamic sampling algorithms for multi-stage  stochastic programs with risk aversion. European Journal of Operational Research  218, 470–483. doi:10.1016/j.ejor.2011.10.056"
},

{
    "location": "readings.html#SDDP.jl-1",
    "page": "Readings",
    "title": "SDDP.jl",
    "category": "section",
    "text": "Two papers about SDDP.jlDowson, O., Kapelevich, L. (2017). SDDP.jl: a Julia package for Stochastic   Dual Dynamic Programming. Optimization Online. link\nDownward, A., Dowson, O., and Baucke, R. (2018). On the convergence of a   cutting plane method for multistage stochastic programming problems with   stagewise dependent price uncertainty. Optimization Online. link"
},

{
    "location": "readings.html#Julia-1",
    "page": "Readings",
    "title": "Julia",
    "category": "section",
    "text": "The organisation\'s websitehttps://julialang.org/The paper describing JuliaBezanson, J., Edelman, A., Karpinski, S., Shah, V.B., 2017. Julia: A Fresh  Approach to Numerical Computing. SIAM Review 59, 65–98. doi:10.1137/141000671"
},

{
    "location": "readings.html#JuMP-1",
    "page": "Readings",
    "title": "JuMP",
    "category": "section",
    "text": "Source code on Githubhttps://www.github.com/JuliaOpt/JuMP.jlThe paper describing JuMPDunning, I., Huchette, J., Lubin, M., 2017. JuMP: A Modeling Language for  Mathematical Optimization. SIAM Review 59, 295–320. doi:10.1137/15M1020575"
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
    "category": "type",
    "text": "SDDPModel(;kwargs...) do ...\n\nend\n\nDescription\n\nThis function constructs an SDDPModel. SDDPModel takes the following keyword arguments. Some are required, and some are optional.\n\nRequired Keyword arguments\n\nstages::Int\n\nThe number of stages in the problem. A stage is defined as each step in time at  which a decion can be made. Defaults to 1.\n\nobjective_bound\n\nA valid bound on the initial value/cost to go. i.e. for maximisation this may  be some large positive number, for minimisation this may be some large negative  number. Users can pass either a single value (which bounds the cost-to-go in all  stages), or a vector of values (one for each stage), or a vector (one element  for each stage) of vectors of values (one value for each markov state in the stage).\n\nsolver::MathProgBase.AbstractMathProgSolver\n\nMathProgBase compliant solver that returns duals from a linear program. If this  isn\'t specified then you must use JuMP.setsolver(sp, solver) in the stage  definition.\n\nOptional Keyword arguments\n\nsense\n\nMust be either :Max or :Min. Defaults to :Min.\n\ncut_oracle::SDDP.AbstractCutOracle\n\nThe cut oracle is responsible for collecting and storing the cuts that define  a value function. The cut oracle may decide that only a subset of the total  discovered cuts are relevant, which improves solution speed by reducing the  size of the subproblems that need solving. Currently must be one of     * DefaultCutOracle() (see DefaultCutOracle for explanation)     * LevelOneCutOracle()(see LevelOneCutOracle for explanation)\n\nrisk_measure\n\nIf a single risk measure is given (i.e. risk_measure = Expectation()), then  this measure will be applied to every stage in the problem. Another option is  to provide a vector of risk measures. There must be one element for every  stage. For example:\n\nrisk_measure = [ NestedAVaR(lambda=0.5, beta=0.25), Expectation() ]\n\nwill apply the i\'th element of risk_measure to every Markov state in the i\'th stage. The last option is to provide a vector (one element for each stage) of vectors of risk measures (one for each Markov state in the stage). For example:\n\nrisk_measure = [\n# Stage 1 Markov 1 # Stage 1 Markov 2 #\n    [ Expectation(), Expectation() ],\n    # ------- Stage 2 Markov 1 ------- ## ------- Stage 2 Markov 2 ------- #\n    [ NestedAVaR(lambda=0.5, beta=0.25), NestedAVaR(lambda=0.25, beta=0.3) ]\n    ]\n\nNote that even though the last stage does not have a future cost function associated with it (as it has no children), we still have to specify a risk measure. This is necessary to simplify the implementation of the algorithm.\n\nFor more help see NestedAVaR or Expectation.\n\nmarkov_transition\n\nDefine the transition probabilties of the stage graph. If a single array is given, it is assumed that there is an equal number of Markov states in each stage and the transition probabilities are stage invariant. Row indices represent the Markov state in the previous stage. Column indices represent the Markov state in the current stage. Therefore:\n\nmarkov_transition = [0.1 0.9; 0.8 0.2]\n\nis the transition matrix when there is 10% chance of transitioning from Markov state 1 to Markov state 1, a 90% chance of transitioning from Markov state 1 to Markov state 2, an 80% chance of transitioning from Markov state 2 to Markov state 1, and a 20% chance of transitioning from Markov state 2 to Markov state 2.\n\nReturns\n\nm: the SDDPModel\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.@state",
    "page": "Reference",
    "title": "SDDP.@state",
    "category": "macro",
    "text": "@state(sp, stateleaving, stateentering)\n\nDescription\n\nDefine a new state variable in the subproblem sp.\n\nArguments\n\nsp               the subproblem\nstateleaving     any valid JuMP @variable syntax to define the value of the state variable at the end of the stage\nstateentering    any valid JuMP @variable syntax to define the value of the state variable at the beginning of the stage\n\nExamples\n\n@state(sp, 0 <= x[i=1:3] <= 1, x0==rand(3)[i] )\n@state(sp,      y        <= 1, y0==0.5        )\n@state(sp,      z            , z0==0.5        )\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.@states",
    "page": "Reference",
    "title": "SDDP.@states",
    "category": "macro",
    "text": "@states(sp, begin\n    stateleaving1, stateentering1\n    stateleaving2, stateentering2\nend)\n\nDescription\n\nDefine a new state variables in the subproblem sp.\n\nArguments\n\nsp               the subproblem\nstateleaving     any valid JuMP @variable syntax to define the value of the state variable at the end of the stage\nstateentering    any valid JuMP @variable syntax to define the value of the state variable at the beginning of the stage\n\nUsage\n\n@states(sp, begin\n    0 <= x[i=1:3] <= 1, x0==rand(3)[i]\n         y        <= 1, y0==0.5\n         z            , z0==0.5\n end)\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.@rhsnoise",
    "page": "Reference",
    "title": "SDDP.@rhsnoise",
    "category": "macro",
    "text": "@rhsnoise(sp, rhs, constraint)\n\nDescription\n\nAdd a constraint with a noise in the RHS vector to the subproblem sp.\n\nArguments\n\nsp         the subproblem\nrhs        keyword argument key=value where value is a one-dimensional array containing the noise realisations\nconstraint any valid JuMP @constraint syntax that includes the keyword defined by rhs\n\nExamples\n\n@rhsnoise(sp, i=1:2, x + y <= i )\n@rhsnoise(sp, i=1:2, x + y <= 3 * rand(2)[i] )\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.@rhsnoises",
    "page": "Reference",
    "title": "SDDP.@rhsnoises",
    "category": "macro",
    "text": "@rhsnoises(sp, rhs, begin\n    constraint\nend)\n\nDescription\n\nThe plural form of @rhsnoise similar to the JuMP macro @constraints.\n\nArguments\n\nSee @rhsnoise.\n\nExamples\n\n@rhsnoises(sp, i=1:2, begin\n    x + y <= i\n    x + y <= 3 * rand(2)[i]\nend)\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.setnoiseprobability!",
    "page": "Reference",
    "title": "SDDP.setnoiseprobability!",
    "category": "function",
    "text": "setnoiseprobability!(sp::JuMP.Model, distribution::Vector{Float64})\n\nDescription\n\nSet the probability distribution of the stagewise independent noise in the sp subproblem.\n\nArguments\n\nsp            the subproblem\ndistribution vector containing the probability of each outcome occuring.   Should sum to 1. Defaults to the uniform distribution.\n\nExamples\n\nIf there are two realizations:\n\nsetnoiseprobability!(sp, [0.3, 0.7])\nsetnoiseprobability!(sp, [0.5, 0.6]) will error!\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.@stageobjective",
    "page": "Reference",
    "title": "SDDP.@stageobjective",
    "category": "macro",
    "text": "@stageobjective(sp, kw=noises, objective)\n\nDescription\n\nDefine an objective that depends on the realization of the stagewise noise. objective can be any valid third argument to the JuMP @objective macro (i.e. @objective(sp, Min, objective)) that utilises the variable kw that takes the realizations defined in noises.\n\nExamples\n\n@stageobjective(sp, w=1:2, w * x)\n@stageobjective(sp, i=1:2, w[i]^2 * x)\n@stageobjective(sp, i=1:2, x[i])\n\n\n\n@stageobjective(sp, objective)\n\nDescription\n\nDefine a deterministic objective.\n\nExamples\n\n@stageobjective(sp, x + y)\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.addconstraintnoise!",
    "page": "Reference",
    "title": "SDDP.addconstraintnoise!",
    "category": "function",
    "text": "addconstraintnoise!(sp::JuMP.Model, x::JuMP.Variable, realizations::Vector{T}) where T\n\nAdd a stagewise-independent coefficient term into the model.\n\nImportant:\n\nRequires a solver that has implemented the MathProgBase.changecoeffs! method.\nCannot be used if JuMP.solvehook is set (e.g., by SDDiP.jl).\n\nIf multiple terms have random coefficents, then the length of realizations must be the same in every case. This can also be used with stagewise-independent right-hand side terms (via @rhsnoise), again with equal numbers of realizations.\n\nExample\n\nIf w ∈ [1,2,3] with equal probability, and we want to add the constraint:     @constraint(sp, 2x + w*x <= 1) Then     wx = addconstraintnoise!(sp, x, [1,2,3])     @constraint(m, 2x + wx <= 1)\n\nNote: This is a bit of a hack. Given a term w*x in the model, for example, in a constraint such as 2x + w*x <= 1, we introduce a dummy variable and a dummy constraint so that: dummyconstraint: dummy_variable == w * x. This allows us to deterministically work around the fact that JuMP cannot modify the coefficient matrix.\n\n\n\n"
},

{
    "location": "apireference.html#Communicating-the-problem-to-the-solver-1",
    "page": "Reference",
    "title": "Communicating the problem to the solver",
    "category": "section",
    "text": "SDDPModel\n@state\n@states\n@rhsnoise\n@rhsnoises\nsetnoiseprobability!\n@stageobjective\naddconstraintnoise!"
},

{
    "location": "apireference.html#SDDP.AbstractRiskMeasure",
    "page": "Reference",
    "title": "SDDP.AbstractRiskMeasure",
    "category": "type",
    "text": "AbstractRiskMeasure\n\nDescription\n\nAbstract type for all risk measures.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.modifyprobability!",
    "page": "Reference",
    "title": "SDDP.modifyprobability!",
    "category": "function",
    "text": "modifyprobability!(measure::AbstractRiskMeasure,\n        riskadjusted_distribution,\n        original_distribution::Vector{Float64},\n        observations::Vector{Float64},\n        m::SDDPModel,\n        sp::JuMP.Model\n)\n\nDescription\n\nCalculate the risk-adjusted probability of each scenario using the \'change-of-probabilities\' approach of Philpott, de Matos, and Finardi,(2013). On solving multistage stochastic programs with coherent risk measures. Operations Research 61(4), 957-970.\n\nArguments\n\nmeasure::AbstractRiskMeasure\n\nThe risk measure\n\nriskadjusted_distribution\n\nA new probability distribution\n\noriginal_distribution::Vector{Float64}\n\nThe original probability distribution.\n\nobservations::Vector{Float64}\n\nThe vector of objective values from the next stage  problems (one for each scenario).\n\nm::SDDPModel\n\nThe full SDDP model\n\nsp::JuMP.Model\n\nThe stage problem that the cut will be added to.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.AVaR",
    "page": "Reference",
    "title": "SDDP.AVaR",
    "category": "type",
    "text": "AVaR(beta::Float64)\n\nThe Average Value @ Risk measure. When beta=0, the measure is the is worst-case, when beta=1 the measure is equivalent to expectation.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.ConvexCombination",
    "page": "Reference",
    "title": "SDDP.ConvexCombination",
    "category": "type",
    "text": "ConvexCombination( (weight::Float64, measure::AbstractRiskMeasure) ... )\n\nCreate a weighted combination of risk measures.\n\nExamples\n\nConvexCombination(\n    (0.5, Expectation()),\n    (0.5, AVaR(0.25))\n)\n\nConvex combinations can also be constructed by adding weighted risk measures together as follows:\n\n0.5 * Expectation() + 0.5 * AVaR(0.5)\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.EAVaR",
    "page": "Reference",
    "title": "SDDP.EAVaR",
    "category": "function",
    "text": "EAVaR(;lambda=1.0, beta=1.0)\n\nDescription\n\nA risk measure that is a convex combination of Expectation and Average Value @ Risk (also called Conditional Value @ Risk).\n\nλ * E[x] + (1 - λ) * AV@R(1-β)[x]\n\nKeyword Arguments\n\nlambda\n\nConvex weight on the expectation ((1-lambda) weight is put on the AV@R component. Inreasing values of lambda are less risk averse (more weight on expecattion)\n\nbeta\n\nThe quantile at which to calculate the Average Value @ Risk. Increasing values  of beta are less risk averse. If beta=0, then the AV@R component is the  worst case risk measure.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.Expectation",
    "page": "Reference",
    "title": "SDDP.Expectation",
    "category": "type",
    "text": "Expectation()\n\nThe expectation risk measure.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.DRO",
    "page": "Reference",
    "title": "SDDP.DRO",
    "category": "type",
    "text": "DRO(radius::Float64)\n\nThe distributionally robust SDDP risk measure. Constructs a DRO risk measure object that allows probabilities to deviate by radius away from the uniform distribution.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.WorstCase",
    "page": "Reference",
    "title": "SDDP.WorstCase",
    "category": "type",
    "text": "WorstCase()\n\nThe worst-case risk measure.\n\n\n\n"
},

{
    "location": "apireference.html#Risk-Measures-1",
    "page": "Reference",
    "title": "Risk Measures",
    "category": "section",
    "text": "AbstractRiskMeasure\nmodifyprobability!\nAVaR\nConvexCombination\nEAVaR\nExpectation\nDRO\nWorstCase"
},

{
    "location": "apireference.html#SDDP.AbstractCutOracle",
    "page": "Reference",
    "title": "SDDP.AbstractCutOracle",
    "category": "type",
    "text": "AbstractCutOracle\n\nDescription\n\nAbstract type for all cut oracles.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.storecut!",
    "page": "Reference",
    "title": "SDDP.storecut!",
    "category": "function",
    "text": "storecut!(oracle::AbstactCutOracle, m::SDDPModel, sp::JuMP.Model, cut::Cut)\n\nDescription\n\nStore the cut cut in the Cut Oracle oracle. oracle will belong to the subproblem sp in the SDDPModel m.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.validcuts",
    "page": "Reference",
    "title": "SDDP.validcuts",
    "category": "function",
    "text": "validcuts(oracle::AbstactCutOracle)\n\nDescription\n\nReturn an iterable list of all the valid cuts contained within oracle.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.allcuts",
    "page": "Reference",
    "title": "SDDP.allcuts",
    "category": "function",
    "text": "allcuts(oracle::AbstactCutOracle)\n\nDescription\n\nReturn an iterable list of all the cuts contained within oracle, not just those that are returned by validcuts.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.DefaultCutOracle",
    "page": "Reference",
    "title": "SDDP.DefaultCutOracle",
    "category": "type",
    "text": "DefaultCutOracle()\n\nDescription\n\nInitialize the default cut oracle.\n\nThis oracle keeps every cut discovered and does not perform cut selection.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.LevelOneCutOracle",
    "page": "Reference",
    "title": "SDDP.LevelOneCutOracle",
    "category": "type",
    "text": "LevelOneCutOracle()\n\nDescription\n\nInitialize the cut oracle for Level One cut selection. See:\n\nV. de Matos, A. Philpott, E. Finardi, Improving the performance of Stochastic Dual Dynamic Programming, Journal of Computational and Applied Mathematics 290 (2015) 196–208.\n\n\n\n"
},

{
    "location": "apireference.html#Cut-Oracles-1",
    "page": "Reference",
    "title": "Cut Oracles",
    "category": "section",
    "text": "AbstractCutOracle\nstorecut!\nvalidcuts\nallcuts\nDefaultCutOracle\nLevelOneCutOracle"
},

{
    "location": "apireference.html#SDDP.StaticPriceInterpolation",
    "page": "Reference",
    "title": "SDDP.StaticPriceInterpolation",
    "category": "type",
    "text": "StaticPriceInterpolation(; kwargs...)\n\nConstuctor for the static price interpolation value function described in\n\nGjelsvik, A., Belsnes, M., and Haugstad, A., (1999). An Algorithm for Stochastic Medium Term Hydro Thermal Scheduling under Spot Price Uncertainty. In PSCC: 13th Power Systems Computation Conference : Proceedings P. 1328. Trondheim: Executive Board of the 13th Power Systems Computation Conference, 1999.\n\nKeyword arguments\n\ndynamics: a function that takes two arguments      1. price: a Float64 that gives the price in the previous stage.      2. noise: a single NoiseRealization of the price noise observed at          the start of the stage.      The function should return a Float64 of the price for the current stage.\ninitial_price: a Float64 for the an initial value for each dimension of the price states.\nrib_locations: an AbstractVector{Float64} giving the points at which to  discretize the price dimension.\nnoise: a finite-discrete distribution generated by DiscreteDistribution\ncut_oracle: any AbstractCutOracle\n\nExample\n\nStaticPriceInterpolation(\n    dynamics = (price, noise) -> begin\n            return price + noise - t\n        end,\n    initial_price = 50.0\n    rib_locations = 0.0:10.0:100.0,\n    noise = DiscreteDistribution([-10.0, 40.0], [0.8, 0.2]),\n)\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.DynamicPriceInterpolation",
    "page": "Reference",
    "title": "SDDP.DynamicPriceInterpolation",
    "category": "type",
    "text": "DynamicPriceInterpolation(; kwargs...)\n\nConstuctor for the dynamic price interpolation value function described in Downward, A., Dowson, O., and Baucke, R. (2018). On the convergence of a cutting plane method for multistage stochastic programming problems with stagewise dependent price uncertainty. Optimization Online.\n\nKeyword arguments\n\ndynamics: a function that takes two arguments      1. price: a Float64 (if uni-variate) or NTuple{N,Float64} if multi-variate          that gives the price in the previous stage.      2. noise: a single NoiseRealization of the price noise observed at          the start of the stage.      The function should return a Float64 (if uni-variate) or NTuple{N,Float64}      if multi-variate that gives the price for the current stage.\ninitial_price: a Float64 (if uni-variate) or NTuple{N,Float64} if multi-variate  that gives the an initial value for each dimension of the price states.\nmin_price: a Float64 (if uni-variate) or NTuple{N,Float64} if multi-variate  that gives the minimum value of each dimension of the price states.\nmax_price: a Float64 (if uni-variate) or NTuple{N,Float64} if multi-variate  that gives the maximum value of each dimension of the price states.\nnoise: a finite-discrete distribution generated by DiscreteDistribution\nlipschitz_constant: the maximum absolute subgradient of any price dimension  in the domain bounded by min_price and max_price.\ncut_oracle: a DynamicCutOracle.\n\nExamples\n\nA uni-variate price process:\n\nDynamicPriceInterpolation(\n    dynamics = (price, noise) -> begin\n            return price + noise - t\n        end,\n    initial_price = 50.0\n    min_price = 0.0,\n    max_price = 100.0,\n    noise = DiscreteDistribution([-10.0, 40.0], [0.8, 0.2]),\n    lipschitz_constant = 1e4\n)\n\nA multi-variate price process:\n\nDynamicPriceInterpolation(\n    dynamics = (price, noise) -> begin\n            return (noise * price[1], noise * price[2] - noise)\n        end,\n    initial_price = (50.0, 50.0),\n    min_price = (0.0,0.0),\n    max_price = (100.0,100.0),\n    noise = DiscreteDistribution([0.5, 1.2], [0.8, 0.2]),\n    lipschitz_constant = 1e4\n)\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.DiscreteDistribution",
    "page": "Reference",
    "title": "SDDP.DiscreteDistribution",
    "category": "type",
    "text": "DiscreteDistribution{T}(observations::AbstractVector{T}, probabilities::AbstractVector{Float64})\n\nCreate a finite discrete distribution of observations supported with probability probabilities.\n\n\n\nDiscreteDistribution{T}(observations::AbstractVector{T})\n\nCreate a finite discrete distribution of observations supported with uniform probability.\n\n\n\n"
},

{
    "location": "apireference.html#Price-Interpolation-1",
    "page": "Reference",
    "title": "Price Interpolation",
    "category": "section",
    "text": "StaticPriceInterpolation\nDynamicPriceInterpolation\nDiscreteDistribution"
},

{
    "location": "apireference.html#JuMP.solve",
    "page": "Reference",
    "title": "JuMP.solve",
    "category": "function",
    "text": "solve(m::SDDPModel; kwargs...)\n\nDescription\n\nSolve the SDDPModel m using SDDP. Accepts a number of keyword arguments to control the solution process.\n\nPositional arguments\n\nm: the SDDPModel to solve\n\nKeyword arguments\n\niteration_limit::Int:  The maximum number of cuts to add to a single stage problem before terminating.\ntime_limit::Real:  The maximum number of seconds to compute for before termination.  Defaults to Inf.\nsimulation::MonteCarloSimulation: see MonteCarloSimulation\nbound_stalling::BoundStalling: see BoundStalling\ncut_selection_frequency::Int:  Frequency (by iteration) with which to rebuild subproblems using a subset of  cuts. Frequent cut selection (i.e. cut_selection_frequency is small) reduces  the size of the subproblems that are solved, but incurrs the overhead of rebuilding  the subproblems. However, infrequent cut selection (i.e.  cut_selection_frequency is large) allows the subproblems to grow large (many  constraints) leading to an increase in the solve time of individual subproblems.  Defaults to 0 (never run).\nprint_level::Int:   0 - off: nothing logged.   1 - on: solve iterations logged.   2 - verbose: detailed timing information is also logged.   Defaults to 1\nlog_file::String:  Relative filename to write the log to disk. Defaults to \"\" (no log written)\nsolve_type:  One of\nAsynchronous() - solve using a parallelised algorithm\nSerial() - solve using a serial algorithm\nDefault chooses automatically based on the number of available processors.\nreduce_memory_footprint::Bool:  Implements the idea proposed in https://github.com/JuliaOpt/JuMP.jl/issues/969#issuecomment-282191105  to reduce the memory consumption when running SDDP. This is an issue if you  wish to save the model m to disk since it discards important information.  Defaults to false.\ncut_output_file::String:  Relative filename to write discovered cuts to disk. Defaults to \"\" (no cuts written)\n\nReturns\n\nstatus::Symbol:  Reason for termination. One of\n:solving\n:interrupted\n:converged\n:iteration_limit\n:bound_stalling\n:time_limit\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.MonteCarloSimulation",
    "page": "Reference",
    "title": "SDDP.MonteCarloSimulation",
    "category": "type",
    "text": "MonteCarloSimulation(;kwargs...)\n\nDescription\n\nCollection of settings to control the simulation phase of the SDDP solution process.\n\nArguments\n\nfrequency::Int\n\nThe frequency (by iteration) with which to run the policy simulation phase of the algorithm in order to construct a statistical bound for the policy. Defaults to 0 (never run).\n\nmin::Float64\n\nMinimum number of simulations to conduct before constructing a confidence interval for the bound. Defaults to 20.\n\nstep::Float64\n\nNumber of additional simulations to conduct before constructing a new confidence interval for the bound. Defaults to 1.\n\nmax::Float64\n\nMaximum number of simulations to conduct in the policy simulation phase. Defaults to min.\n\nconfidence::Float64\n\nConfidence level of the confidence interval. Defaults to 0.95 (95% CI).\n\ntermination::Bool\n\nWhether to terminate the solution algorithm with the status :converged if the deterministic bound is with in the statistical bound after max simulations. Defaults to false.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.BoundStalling",
    "page": "Reference",
    "title": "SDDP.BoundStalling",
    "category": "type",
    "text": "BoundStalling(;kwargs...)\n\nDescription\n\nCollection of settings to control the bound stalling convergence test.\n\nArguments\n\niterations::Int\n\nTerminate if the maximum deviation in the deterministic bound from the mean over the last iterations number of iterations is less than rtol (in relative terms) or atol (in absolute terms).\n\nrtol::Float64\n\nMaximum allowed relative deviation from the mean. Defaults to 0.0\n\natol::Float64\n\nMaximum allowed absolute deviation from the mean. Defaults to 0.0\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.Asynchronous",
    "page": "Reference",
    "title": "SDDP.Asynchronous",
    "category": "type",
    "text": "Asynchronous(; kwargs...)\n\nDefine\n\nType used to dispatch and control the behaviour of the asynchronous solution algorithm.\n\nArguments\n\nslaves::Vector{Int} the pid\'s of the slave processes. Defaults to workers()\nstep::Float64 the number of iterations to complete before adding another  slave. Used to replicated the scenario incrementation behaviour of  V. de Matos,A. Philpott, E. Finardi, Improving the performance of Stochastic  Dual Dynamic Programming, Journal of Computational and Applied Mathematics  290 (2015) 196–208.\n\nExamples\n\nAsynchronous() # load on all workers\nAsynchronous(slaves=[2,3,4]) # load slaves on processes 2, 3, and 4\nAsynchronous(step=10) # perform 10 iterations before adding new slave\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.Serial",
    "page": "Reference",
    "title": "SDDP.Serial",
    "category": "type",
    "text": "Serial()\n\nDefine\n\nType used to dispatch the serial solution algorithm\n\n\n\n"
},

{
    "location": "apireference.html#Solving-the-problem-efficiently-1",
    "page": "Reference",
    "title": "Solving the problem efficiently",
    "category": "section",
    "text": "solve\nMonteCarloSimulation\nBoundStalling\nAsynchronous\nSerial"
},

{
    "location": "apireference.html#SDDP.simulate",
    "page": "Reference",
    "title": "SDDP.simulate",
    "category": "function",
    "text": "simulate(m::SDDPPModel,variables::Vector{Symbol};\n    noises::Vector{Int}, markovstates::Vector{Int})\n\nDescription\n\nPerform a historical simulation of the current policy in model  m.\n\nnoises is a vector with one element for each stage giving the index of the (in-sample) stagewise independent noise to sample in each stage. markovstates is a vector with one element for each stage giving the index of the (in-sample) markov state to sample in each stage.\n\nExamples\n\nsimulate(m, [:x, :u], noises=[1,2,2], markovstates=[1,1,2])\n\n\n\nresults = simulate(m::SDDPPModel, N::Int, variables::Vector{Symbol})\n\nDescription\n\nPerform N Monte-Carlo simulations of the current policy in model m saving the values of the variables named in variables at every stage.\n\nresults is a vector containing a dictionary for each simulation. In addition to the variables specified in the function call, other special keys are:\n\n:stageobjective - costs incurred during the stage (not future)\n:obj            - objective of the stage including future cost\n:markov         - index of markov state visited\n:noise          - index of noise visited\n:objective      - Total objective of simulation\n\nAll values can be accessed as follows\n\nresults[simulation index][key][stage]\n\nwith the exception of :objective which is just\n\nresults[simulation index][:objective]\n\nExamples\n\nresults = simulate(m, 10, [:x, :u])\nresults[1][:objective] # objective of simulation 1\nmean(r[:objective] for r in results) # mean objective of the simulations\nresults[2][:x][3] # value of :x in stage 3 in second simulation\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.getbound",
    "page": "Reference",
    "title": "SDDP.getbound",
    "category": "function",
    "text": "getbound(m)\n\nDescription\n\nGet the lower (if minimizing), or upper (if maximizing) bound of the solved SDDP model m.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.newplot",
    "page": "Reference",
    "title": "SDDP.newplot",
    "category": "function",
    "text": "SDDP.newplot()\n\nDescription\n\nInitialize a new SimulationPlot.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.addplot!",
    "page": "Reference",
    "title": "SDDP.addplot!",
    "category": "function",
    "text": "SDDP.addplot!(p::SimulationPlot, ivals::AbstractVector{Int}, tvals::AbstractVector{Int}, f::Function; kwargs...)\n\nDescription\n\nAdd a new figure to the SimulationPlot p, where the y-value is given by f(i, t) for all i in ivals (one for each series) and t in tvals (one for each stage).\n\nKeywords\n\nxlabel: set the xaxis label;\nylabel: set the yaxis label;\ntitle: set the title of the plot;\nymin: set the minimum y value;\nymax: set the maximum y value;\ncumulative: plot the additive accumulation of the value across the stages;\ninterpolate: interpolation method for lines between stages. Defaults to \"linear\"  see the d3 docs\n\nfor all options.\n\nExamples\n\nresults = simulate(m, 10)\np = SDDP.newplot()\nSDDP.addplot!(p, 1:10, 1:3, (i,t)->results[i][:stageobjective][t])\n\n\n\n"
},

{
    "location": "apireference.html#Base.show-Tuple{SDDP.SimulationPlot}",
    "page": "Reference",
    "title": "Base.show",
    "category": "method",
    "text": "show(p::SimulationPlot)\n\nLaunch a browser and render the SimulationPlot plot p.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.plotvaluefunction",
    "page": "Reference",
    "title": "SDDP.plotvaluefunction",
    "category": "function",
    "text": " SDDP.plotvaluefunction(m::SDDPModel, stage::Int, markovstate::Int, states::Union{Float64, AbstractVector{Float64}}...; label1=\"State 1\", label2=\"State 2\")\n\nDescription\n\nPlot the value function of stage stage and Markov state markovstate in the SDDPModel m at the points in the discretized state space given by states. If the value in states is a real number, the state is evaluated at that point. If the value is a vector, the state is evaluated at all the points in the vector. At most two states can be vectors.\n\nExamples\n\nSDDP.plotvaluefunction(m, 2, 1, 0.0:0.1:1.0, 0.5, 0.0:0.1:1.0; label1=\"State 1\", label2=\"State 3\")\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.getsubproblem",
    "page": "Reference",
    "title": "SDDP.getsubproblem",
    "category": "function",
    "text": "getsubproblem(m::SDDPModel, t::Int, i::Int)\n\nGet the subproblem in stage t and Markov state i from the SDDPModel m.\n\n\n\n"
},

{
    "location": "apireference.html#Understanding-the-solution-1",
    "page": "Reference",
    "title": "Understanding the solution",
    "category": "section",
    "text": "simulate\ngetbound\nnewplot\naddplot!\nshow(::SDDP.SimulationPlot)\nplotvaluefunction\ngetsubproblem"
},

{
    "location": "apireference.html#SDDP.writecuts!",
    "page": "Reference",
    "title": "SDDP.writecuts!",
    "category": "function",
    "text": "writecuts!(filename::String, m::SDDPModel; onlyvalid=false)\n\nWrites all cuts from model m to filename.\n\nIf onlyvalid is true, write the cuts returned from validcuts, else write the cuts returned from allcuts.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.loadcuts!",
    "page": "Reference",
    "title": "SDDP.loadcuts!",
    "category": "function",
    "text": "loadcuts!(m::SDDPModel, filename::String)\n\nLoad cuts from the file created using the cut_output_file argument in solve.\n\nExample\n\nm = SDDPModel() do ... end\nstatus = solve(m; cut_output_file=\"path/to/m.cuts\")`\nm2 = SDDPModel() do ... end\nloadcuts!(m2, \"path/to/m.cuts\")\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.savemodel!",
    "page": "Reference",
    "title": "SDDP.savemodel!",
    "category": "function",
    "text": "SDDP.savemodel!(filename::String, m::SDDPModel)\n\nSave the SDDPModel m to the location filename. Can be loaded at a later date with m = SDDP.loadmodel(filename).\n\nNote: this function relies in the internal Julia Base.serializefunction. It should not be relied on to save an load models between versions of Julia (i.e between v0.5 and v0.6). For a longer-term solution, see SDDP.loadcuts! for help.\n\n\n\n"
},

{
    "location": "apireference.html#SDDP.loadmodel",
    "page": "Reference",
    "title": "SDDP.loadmodel",
    "category": "function",
    "text": "loadmodel(filename::String)\n\nLoad a model from the location filename that was saved using SDDP.savemodel!.\n\nNote: this function relies in the internal Julia Base.serializefunction. It should not be relied on to save an load models between versions of Julia (i.e between v0.5 and v0.6). For a longer-term solution, see SDDP.loadcuts! for help.\n\n\n\n"
},

{
    "location": "apireference.html#Read-and-write-the-model-to-disk-1",
    "page": "Reference",
    "title": "Read and write the model to disk",
    "category": "section",
    "text": "writecuts!\nloadcuts!\nsavemodel!\nloadmodel"
},

]}
