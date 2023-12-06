using JuMP
import Ipopt
import Plots

function solve_v_inner(f, x̅, λ)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, -2 <= x <= 2, start = x̅)
    @objective(model, Min, f(x) - λ * (x - x̅))
    optimize!(model)
    return objective_value(model), value(x)
end

function solve_extreme_z(f, x̅, λ)
    ret = [solve_v_inner(fi, x̅, λ) for fi in f]
    _, i = findmin(first, ret)
    return ret[i]
end

function solve_average_z(f, x̅, λ)
    f_average(x) = sum(fi(x) for fi in f) / length(f)
    return solve_v_inner(f_average, x̅, λ)
end

function solve_outer(
    f::Vector{<:Function},
    x̅;
    plot::Bool = false,
    atol::Float64 = 1e-5,
    iteration_limit::Int = 40,
    solve_inner::Function = solve_extreme_z,
)
    upper_bound = minimum(fi(x̅) for fi in f)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, λ, start = 0)
    @variable(model, t <= upper_bound)
    @objective(model, Max, t)
    if plot
        plot_x = Plots.plot(xlims = (-1.5, 1.5), ylims = (0, 10), xlabel = "x")
        plot_λ = Plots.plot(ylims = (0, upper_bound), xlabel = "λ")
        Plots.vline!(plot_x, [x̅]; color = :orange, label = "x")
    end
    X = -2:0.01:2
    if plot
        for fi in f
            Plots.plot!(plot_x, fi, X; label = false)
        end
    end
    lb, ub, old_l, k = -Inf, Inf, 0.0, 0
    while k < iteration_limit && abs(ub - lb) > atol
        optimize!(model)
        l = value(λ)
        if abs(old_l - l) > 1
            l = (old_l + l) / 2  # Some minor regularizationn
        end
        old_l = l
        L, x = solve_inner(f, x̅, l)
        lb, ub = max(L, lb), value(t)
        println("l = $l  lb = $(round(lb; digits = 3))  ub = $(round(ub; digits = 3))")
        @constraint(model, t <= L + (x̅ - x) * (λ - l))
        if plot
            Plots.plot!(plot_x, x -> L + (x - x̅) * l, X;  label = false, color = :gray)
            Plots.plot!(plot_λ, λ -> L + (x̅ - x) * (λ - l), X;  label = false, color = :gray)
        end
        k += 1
    end
    optimize!(model)
    if plot
        Plots.vline!(plot_λ, [value(λ)]; color = :orange, label = "λ")
        display(Plots.plot(plot_x, plot_λ))
    end
    return value(λ)
end

f = [x -> 1/10 + (x - 1)^2, x -> 1/5 + (x + 1)^2]
λ = solve_outer(f, 0.0; plot = true)
λ = solve_outer(f, 0.0; plot = true, solve_inner = solve_average_z)
