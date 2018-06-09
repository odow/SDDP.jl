#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    Example: newsvendor.

    This example is based on the classical newsvendor problem, but features
    an AR(1) spot-price.

    V(x[t-1], ω[t]) =         max p[t] × u[t]
                       subject to x[t] = x[t-1] - u[t] + ω[t]
                                  u[t] ∈ [0, 1]
                                  x[t] ≥ 0
                                  log(p[t]) = log(p[t-1]) + log(ϕ[t])

    x[0] = 2.0
    p[0] = 1.5
    ω[t] ~ {0, 0.05, 0.10, ..., 0.45, 0.5} with uniform probability.
    ϕ[t] ~ {0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1} with uniform probability.
=#

using SDDP, JuMP, Clp, Base.Test

function newsvendor_example(DISCRETIZATION = 1)

    srand(10)
    function buildvaluefunction(stage, markovstate)
        ϕ = log.([0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1])
        INITIAL_PRICE = log(1.50)
        MIN_PRICE     = INITIAL_PRICE + stage * minimum(ϕ)
        MAX_PRICE     = INITIAL_PRICE + stage * maximum(ϕ)
        NOISES        = DiscreteDistribution(ϕ)

        function pricedynamics(price, noise)
            price + noise
        end

        if DISCRETIZATION == 1
            return DynamicPriceInterpolation(
                dynamics       = pricedynamics,
                initial_price  = INITIAL_PRICE,
                min_price      = MIN_PRICE,
                max_price      = MAX_PRICE,
                noise          = NOISES,
            lipschitz_constant = 10.0
            )
        else
            return StaticPriceInterpolation(
                dynamics       = pricedynamics,
                initial_price  = INITIAL_PRICE,
                rib_locations  =  collect(linspace(MIN_PRICE, MAX_PRICE, DISCRETIZATION)),
                noise          = NOISES
            )
        end
    end

    m = SDDPModel(
        sense             = :Min,
        stages            = 5,
        objective_bound   = -10,
        solver            = ClpSolver(),
        value_function    = buildvaluefunction
                                            ) do sp, t
        @state(sp, x >= 0, x0 == 2)
        @variable(sp, 0 <= u <= 1)
        @rhsnoise(sp, ω = 0:0.05:0.5, x  == x0 - u + ω)
        @stageobjective(sp, price -> -exp(price) * u)
    end
    return m
end

if length(ARGS) > 0
    if ARGS[1] == "static"
        m = newsvendor_example(15)
        srand(123)
        status = SDDP.solve(m,
            iteration_limit = 10
        )
        s = simulate(m, 1_000, [:x, :x0, :u, :price])

        using Plots
        gr()
        Ω = 0:0.05:0.5
        plot(
            SDDP.publicationplot(s, :x, title="Outgoing State", ylabel="Widgets"),
            SDDP.publicationplot(s, x->exp.(x[:price]), title="Price", ylabel="\$/Widget"),
            SDDP.publicationplot(s, x->Ω[x[:noise]], title="Noise", ylabel="Widgets"),
            SDDP.publicationplot(s, :u, title="Control", ylabel="Widgets"),
            layout        = (2,2),
            size          = (1000, 600),
            titlefont     = Plots.font("times", 14),
            guidefont     = Plots.font("times", 14),
            tickfont      = Plots.font("times", 14),
            bottom_margin = 7.5Plots.mm,
            left_margin   = 5Plots.mm
        )
        savefig("static.pdf")

        bounds = Float64[]
        for N in 2:20
            m = newsvendor_example(N)
            srand(123)
            status = SDDP.solve(m,
                iteration_limit = 50,
                print_level=0
            )
            push!(bounds , getbound(m))
        end
        open("static.bound", "w") do io
            for b in bound
                println(io, b)
            end
        end

    elseif ARGS[1] == "dynamic"
        m = newsvendor_example(1)
        srand(123)
        status = SDDP.solve(m,
            iteration_limit = 500
        )
        bound = [l.bound for l in m.log]
        open("dynamic.bound.csv", "w") do io
            for (i,b) in enumerate(bound)
                println(io, "$(i), $(b)")
            end
        end
    elseif ARGS[1] == "hybrid"
        # warm up Julia JIT
        m = newsvendor_example(1)
        SDDP.solve(m, iteration_limit = 1)
        m = newsvendor_example(3)
        SDDP.solve(m, iteration_limit = 1)

        # Dynamic Only
        m = newsvendor_example(1)
        srand(123)
        status = SDDP.solve(m,
            time_limit=20.0
        )
        sol_time = [l.timetotal for l in m.log]
        bound = [l.bound for l in m.log]
        open("dynamic.1.csv", "w") do io
            for (t, b) in zip(sol_time, bound)
                println(io, "$(t), $(b)")
            end
        end

        # Static Only
        m = newsvendor_example(5)
        srand(123)
        status = SDDP.solve(m,
            time_limit=20.0
        )
        sol_time = [l.timetotal for l in m.log]
        bound = [l.bound for l in m.log]
        open("static.5.csv", "w") do io
            for (t, b) in zip(sol_time, bound)
                println(io, "$(t), $(b)")
            end
        end

        # Hybrid
        m = newsvendor_example(5)
        srand(123)
        status = SDDP.solve(m,
            bound_stalling = BoundStalling(iterations = 5, atol=1e-4),
            cut_output_file = "hybrid.cuts"
        )

        sol_time = [l.timetotal for l in m.log]
        bound = [l.bound for l in m.log]
        offset = maximum(sol_time)
        m = newsvendor_example(1)
        SDDP.loadcuts!(m, "hybrid.cuts")
        rm("hybrid.cuts")
        status = SDDP.solve(m,
            time_limit=20.0 - offset
        )
        for l in m.log
            push!(sol_time, l.timetotal + offset)
            push!(bound, l.bound)
        end
        open("hybrid.5.csv", "w") do io
            for (t, b) in zip(sol_time, bound)
                println(io, "$(t), $(b)")
            end
        end
    end
end
