# Tutorial Eleven: Infinite-Horizon Stochastic Dynamic Programming

This tutorial discusses discusses Infinite-Horizon Stochastic Dynamic Programming. The implementation of infinite-horizon SDDP uses the average cost method. Another possible method could be using the discounted cost method however this method converges much slower than the average cost method. Infinite-horizon SDDP finds the optimal steady state policy of the problem.

My Honors [thesis](https://github.com/shasafoster/SDDP.jl/blob/master/docs/src/assets/foster_thesis.pdf) may be useful for further understanding the underlying theory of Infinite-Horizon Stochastic Dynamic Programming. Ben Fulton applied infinite-horizon SDDP when modelling various scenarios in the New Zealand electricity market thus his [thesis](https://github.com/shasafoster/SDDP.jl/blob/master/docs/src/assets/fulton_thesis.pdf) may also be of interest. 

The differencs between infinite-horizon SDDP and SDDP will be explained by formulating and solving a multistage stochastic dynamic program with each method. 

## The Problem

We will continue using the hydrothermal scheduling problem, the most common application of stochastic dual dynamic programming. 

The difference between this formulation and the formulation used in previous tutorials is the presence of a terminal cost-to-go function. The terminal cost-to-go is an important part of more developed hydrothermal scheduling problems. 

In the context of hydrothermal sceduling the terminal cost is based on the concept of a marginal value of water. Water at the end-of-horizon has a value because it can be used to generate electricity. The water is said to have a *marginal* value because a m^3 of water is worth more to when the reservoir is empty compared to when our reseroivr is full.

In the absence of a terminal cost-to-go function in the context of the hydrothermal sceduling problem the policy would leave the reservoir empty at the end of the final stage. As water in hydro reservoirs has value this outcome is undesirable and would not happen in reality. A terminal cost-to-go function is used to penalise such action. A terminal cost, which is a function of the final reserovir levels at the end of the final stage is included in the terminal stage objective. 

The terminal cost function used in the simple example is shown in the chart below:

![terminal_cost_chart](../assets/terminal_cost_chart.png)

From the chart we see 



![marginal_water_value_chart](../assets/marginal_water_value_chart.png)


## Formulating and solving the problem in SDDP

```julia
using SDDP, JuMP, Clp
m = SDDPModel(
                 sense = :Min,
                 stages = 12,
                 solver = ClpSolver(),
                 objective_bound = 0.0,
                 is_infinite = false,
                                        ) do sp, t
                                        
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        terminalcost       >= 0
     end)
     
    inflow = [50.0, 50.0, 50.0]
    fuel_cost = [50.0, 100.0, 150.0]
    terminal_marginal_cost = [-30 -20 -10]
    intercept = [200 100 50]
    
    @expression(md, terminalcost
        for i in 1:length(marginal_cost)
	    terminalcost >= terminal_marginal_cost[i] * outgoing_volume + intercept[i])
	end
    )
    
    @constraints(sp, begin
        incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume
        thermal_generation + hydro_generation == 150
    end)
    
    if t < 12
        @stageobjective(sp, fuel_cost[t] * thermal_generation)
    elseif t == 12
        # Final stage has extra term in objectivve - the terminal cost
        @stageobjective(sp, fuel_cost[t] * thermal_generation + terminalcost)
    end
end
```

To solve this problem, we use the solve method:

status = solve(m; iteration_limit=5)

## How the formulation changes
In formulating many stochastic dynamic programs, a terminating cost-to-go function is necessary. However, this terminating cost-to-go function is an assumption of many fomulations. Solving a multi-stage stochastic dynamic problem with infinite-horizon stochastic dynamic programming (infinite-horizon SDDP), eliminates the need for a terminating cost-to-go function. 


This concludes our tutorial 12 for SDDP.jl on infinite-horizon SDDP. 
