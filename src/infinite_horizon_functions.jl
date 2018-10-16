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
        y = max((sum(vA(new_cuts[i,4+dim_state:end], new_cuts[:,4:3+dim_state]),2) .+ new_cuts[:,3])...)
        delta = min((y - (sum(vA(new_cuts[i,4+dim_state:end], current_cuts[:,4:3+dim_state]),2) .+ current_cuts[:,3]))...)
        if delta > 0
            append!(delta_arr, delta)
        end
    end

    Delta = min(delta_arr...)
    new_cuts[:,3] -= Delta
    return new_cuts, delta_arr
end



function convergence_test(allcuts_fp::String, ub_arr::Array{Float64,1})
    "Function determines area under 1D approximation of all of stage 1 cuts"
    # Total runtime about 1.1 second

    ub_arr = [242544.0, 84862.0, 82319.0, 42345.0, 150187.0, 137876.0, 5724.0]
    cuts = cutsout[get_idx(cutsout, 1),:]
    dim_state = div(size(cuts,2)-3,2)
    sum_state = sum(cuts[:,4+dim_state:end],2)
    idx = cuts[:,3] .> repeat([1e-6], inner=size(cuts,1))
    weighted_product = sum(cuts[:,4:3+dim_state] .* cuts[:,4+dim_state:end],2)
    betas = weighted_product[idx] ./ sum_state[idx]
    alphas = cuts[idx,3]
    idx = .!(isnan.(alphas) .| isnan.(betas))
    alphas = alphas[idx]
    betas = betas[idx]

    duals = zeros(0)
    dominating_cuts_idx = zeros(0)
    X = collect(linspace(0,sum(ub_arr),10000))

    for x̄ in X
        val, idx = findmax(alphas + betas * x̄)
        append!(dominating_cuts_idx, idx)
    end

    idx = unique(dominating_cuts_idx)
    alphas = alphas[Int.(idx)]
    betas = betas[Int.(idx)]
    alphas = alphas

    # Takes approx 2 seconds
    integral_sum = 0
    x_max = max((-alphas ./ betas)...)
    n_steps = 1000000
    step = x_max/(n_steps+1)
    x = collect(linspace(0,x_max,n_steps))

    alphasM = repmat(alphas, 1, n_steps)
    betasM = repmat(betas, 1, n_steps)
    xM = repmat(x, 1, 4)'
    vals, idx = findmax(alphasM .+ betasM .* xM, 1)

    integral = vals * step * ones(n_steps)
    integral
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
