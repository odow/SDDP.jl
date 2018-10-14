using Gurobi

function write_cuts(filepath::String, cuts::Array{Float64,2}, file_modifier::String)
    open(filepath, file_modifier) do file
        for i in 1:size(cuts,1)
            cut = cuts[i,:]
            write(file, string(Int(cut[1])), ",", string(Int(cut[2])),",", string(cut[3]))
            for v in cut[4:end]
                write(file, ",", string(v))
            end
            write(file,"\n")
        end
    end
end


function save_deltas(delta_fp::String, delta_arr::Array{Float64,1}, file_modifier::String)
    open(delta_fp, file_modifier) do file
        for d in delta_arr
            write(file, string(d), ",")
        end
        write(file,"\n")
    end
end


function shift_newcuts_down!(new_cuts::Array{Float64,2}, stageTcuts_fp::String)
    current_cuts = readcsv(stageTcuts_fp)
    delta_arr = zeros(0)
    dim_state = div(size(new_cuts,2)-3,2)

    # Determine the minimum distance between the new cut surface
    # and the current cut surface for each of the new cuts sampled state
    for i in 1:size(new_cuts,1)
        y = max((sum(vA(new_cuts[i,(4+dim_state):end], new_cuts[:,4:(3+dim_state)]),2) .+ new_cuts[:,3])...)
        delta = min((y - (sum(vA(new_cuts[i,(4+dim_state):end], current_cuts[:,4:(3+dim_state)]),2) .+ current_cuts[:,3]))...)
        if delta > 0
            append!(delta_arr, delta)
        end
    end

    Delta = min(delta_arr...)
    new_cuts[:,3] -= Delta
    return new_cuts, delta_arr
end


function convergence_test(stageTcuts_fp::String)
    "Function determines area under 1D approximation of all of stage 1 cuts"

    theta = [0.000725151, 0.00067598, 0.001131078, 0.000363842, 0.000395587, 0.00026111, 0.000725151]
    scalefactor = 10000
    cuts = readcsv(stageTcuts_fp)
    alphas = zeros(0)
    betas_1D = zeros(0)


    for i in 1:size(cuts,1)
        cut = cuts[i,:]
        if (sum(cut[11:17]) > 0.01) & (cut[3] > 0.01) & (abs(sum(cut[4:10])) > 0.01)
            alpha = cut[3]
            beta = sum((cut[11:17] .* cut[4:10]) ./ theta) / (sum(cut[11:17]) * scalefactor) # Slope
            append!(alphas, alpha)
            append!(betas_1D, beta)
        end
    end

    nb_cuts = length(alphas)
    duals = zeros(0)
    dominating_cuts_indexes = zeros(0)
    x_arr = collect(linspace(1,5000,1000))
    for x_hat in x_arr
        m = Model(solver=GurobiSolver(OutputFlag=0))
        @variable(m, 0 <= theta <= Inf)
        @variable(m, 0 <= x <= Inf)
        @objective(m, Min, 1theta)

        # Reference constraints by cut[i] so we we can see later which cuts are dominating
        @constraint(m, cut[i=1:nb_cuts], betas_1D[i]*x + alphas[i] <= theta)

        # Add constraint to set current energy storage
        # 1000 scalefactor as x_hat is in GWh not MWh
        state_constraint = @constraint(m, x <= 1000 * x_hat)

        status = solve(m)
        if status != :Optimal
            println("LP not optimal for, x:", x)
        else
            # Store marginal value for this energy state
            append!(duals, getdual(state_constraint))

            # If the shadow price is zero the constraint is not binding
            # We only want the dominating cut (which is the binding constraint)
            # for this given x_hat (energy state)
            for i in 1:nb_cuts
                if abs(getdual(cut[i])) > 1e-8
                    append!(dominating_cuts_indexes, Int(i))
                end
            end
        end
    end

    dominating_cuts_indexes = unique(dominating_cuts_indexes)

    alphas = alphas[Int.(dominating_cuts_indexes)] / 1000
    betas_1D = betas_1D[Int.(dominating_cuts_indexes)]

    # Takes approx 2 seconds
    integral_sum = 0
    x_max = max((-alphas ./ betas_1D)...)
    n_steps = 1000000
    for x in collect(linspace(0,x_max,n_steps))
        integral_sum += max(append!(alphas + betas_1D * x, 0)...) * (x_max/(n_steps+1))
    end

    integral_sum
end

function print_covergence_metrics(time_arr::Array{Float64,1},
                                  lb_arr::Array{Float64,1},
                                  Delta_arr::Array{Float64,1},
                                  SDdelta::Array{Float64,1},
                                  terminalcost_integral::Array{Float64,1})
    println("-------------------------------------------------------------------------------")
    println("Iteration  Time (s)   Bound      Delta      sd(delta)    Terminal Cost Integral")
    println(humanize(1), "     ", humanize(time_arr[1],"5.2f"), "  ", humanize(lb_arr[1],"8.3f"))
    for i in 2:length(time_arr)
        println(humanize(i), "     ",
                humanize(time_arr[i],"5.2f"), "  ",
                humanize(lb_arr[i],"8.3f"),"   ",
                humanize(Delta_arr[i-1],"8.3f"), "   ",
                humanize(SDdelta[i-1],"8.3f"), "   ",
                humanize(terminalcost_integral[i-1],"8.3f"))
    end
    println("Net runtime (minutes): ", Int(round(sum(time_arr) / 60,0)))
    println("-------------------------------------------------------------------------------")
end
