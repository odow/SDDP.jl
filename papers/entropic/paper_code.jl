#  Copyright 2020-21, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This code is used to build the figures contained in the paper
#
# On solving multistage stochastic programs with the entropic risk measure.
# Dowson, O., Morton, D.P., and Pagnoncelli, B.K.

using SDDP

import Distributions
import Gurobi
import Random
import Statistics

const GRB_ENV = Gurobi.Env()

function gurobi_optimizer()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    return model
end

function _write_matrix(io, data)
    for i in 1:size(data, 1)
        for j in 1:size(data, 2)
            print(io, data[i, j], " ")
        end
        print(io, "\n")
    end
    return
end

###
### FIGURE 1
###

function _ent(z, p, γ)
    f = [pk * exp(γ * zk) for (pk, zk) in zip(p, z)]
    f ./= sum(f)
    return f
end

function _eoh_cvar(z, p, γ)
    y = zeros(length(z))
    α = 0.0
    for k in sortperm(z; rev = true)
        if α >= (1 - γ)
            break
        end
        y[k] = min(p[k], (1 - γ) - α) / (1 - γ)
        α += y[k] * (1 - γ)
    end
    return y
end

function _nested_cvar(z1, p1, z2, p2, γ)
    y2 = _eoh_cvar(z2, p2, γ)
    y1 = _eoh_cvar(z1 .+ y2' * z2, p1, γ)
    return vcat(y1[1] .* y2, y1[2] .* y2)
end

function build_figure_1()
    @info "Building Figure 1"
    z = collect(1:8)
    p = fill(1 / 8, 8)
    γ = 0.4
    data = hcat(
        z,
        p,
        _ent(z, p, γ),
        _eoh_cvar(z, p, γ),
        _nested_cvar([0, 4], [0.5, 0.5], 1:4, fill(0.25, 4), γ),
    )
    open(joinpath(@__DIR__, "risk_set.dat"), "w") do io
        return _write_matrix(io, data)
    end
    return
end

build_figure_1()

###
### FIGURE 3
###

function _solve_example(risk_measure)
    X = Dict(2.0 => 0.1, 1.8 => 0.9)
    Y = Dict(2.2 => 0.1, 1.7 => 0.9)
    Z = Dict(2.0 => 0.1, 1.0 => 0.9)
    two_stage_model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Min,
        lower_bound = 0.0,
        optimizer = gurobi_optimizer,
    ) do sp, t
        @variable(sp, 0 <= x <= 1, SDDP.State, Bin, initial_value = 0)
        if t == 1
            @stageobjective(sp, 0)
        else
            Ω = [(ω1, ω2, ω3) for (ω1, _) in X for (ω2, _) in Y for (ω3, _) in Z]
            P = [
                p1 * p2 * p3 for (_, p1) in X for (_, p2) in Y for (_, p3) in Z
            ]
            SDDP.parameterize(sp, Ω, P) do ω
                @stageobjective(sp, ω[1] * x.in + ω[2] * (1 - x.in) + ω[3])
            end
        end
    end
    three_stage_model = SDDP.LinearPolicyGraph(
        stages = 3,
        sense = :Min,
        lower_bound = 0.0,
        optimizer = gurobi_optimizer,
    ) do sp, t
        @variable(sp, 0 <= x <= 1, SDDP.State, Bin, initial_value = 0)
        if t == 1
            @stageobjective(sp, 0)
        elseif t == 2
            Ω = [(ω1, ω2) for (ω1, _) in X for (ω2, _) in Y]
            P = [p1 * p2 for (_, p1) in X for (_, p2) in Y]
            SDDP.parameterize(sp, Ω, P) do ω
                @stageobjective(sp, ω[1] * x.in + ω[2] * (1 - x.in))
            end
        elseif t == 3
            SDDP.parameterize(sp, [ω3 for (ω3, _) in Z], [p3 for (_, p3) in Z]) do ω
                @stageobjective(sp, ω)
            end
        end
    end
    SDDP.train(
        two_stage_model,
        iteration_limit = 100,
        print_level = 0,
        risk_measure = risk_measure,
        duality_handler = SDDP.LagrangianDuality(),
    )
    SDDP.train(
        three_stage_model,
        iteration_limit = 100,
        print_level = 0,
        risk_measure = risk_measure,
        duality_handler = SDDP.LagrangianDuality(),
    )
    two_stage_x = SDDP.simulate(two_stage_model, 1, [:x])[1][1][:x].out
    three_stage_x = SDDP.simulate(three_stage_model, 1, [:x])[1][1][:x].out
    return two_stage_x, three_stage_x
end

function build_figure_3()
    @info "Building Figure 3"
    open(joinpath(@__DIR__, "results_avar.dat"), "w") do io
        for γ in 0.5:0.01:1.0
            println("  AVaR(γ)     = $γ")
            two, three = _solve_example(SDDP.AVaR(1 - γ))
            println(io, "$(rpad(γ, 5, '0')) $(Int(two)) $(Int(three))")
        end
    end

    open(joinpath(@__DIR__, "results_entropic.dat"), "w") do io
        for γ in vcat(0.0, 4:0.05:5, 10.0)
            println("  Entropic(γ) = $γ")
            two, three = _solve_example(SDDP.Entropic(γ))
            println(io, "$(lpad(γ, 4)) $(Int(two)) $(Int(three))")
        end
    end
    return
end

build_figure_3()

###
### FIGURE 6
###

mutable struct CyclicHistorical{T,S} <: SDDP.AbstractSamplingScheme
    scenarios::Vector{Vector{Tuple{T,S}}}
    current::Int
end

function SDDP.sample_scenario(
    ::SDDP.PolicyGraph{T},
    sampling_scheme::CyclicHistorical{T,S};
    # Ignore the other kwargs because the user is giving
    # us the full scenario.
    kwargs...,
) where {T,S}
    if sampling_scheme.current > length(sampling_scheme.scenarios)
        sampling_scheme.current = 1
    end
    scenario = sampling_scheme.scenarios[sampling_scheme.current]
    sampling_scheme.current += 1
    return scenario, false
end

function _run_instance(risk_measure)
    Ω, P = [(s = 1.11, b = 1.02), (s = 1.04, b = 1.06)], [0.2, 0.8]
    model = SDDP.LinearPolicyGraph(
        stages = 5,
        sense = :Max,
        upper_bound = 5.0,
        optimizer = gurobi_optimizer,
    ) do sp, t
        set_silent(sp)
        @variable(sp, stocks >= 0, SDDP.State, initial_value = 0)
        @variable(sp, bonds >= 0, SDDP.State, initial_value = 1)
        @variable(sp, consumption >= 0)
        @constraint(
            sp,
            c,
            stocks.out + bonds.out + consumption - stocks.in - bonds.in == 0
        )
        @stageobjective(sp, consumption)
        if t > 1
            SDDP.parameterize(sp, Ω, P) do ω
                set_normalized_coefficient(c, stocks.in, -ω.s)
                return set_normalized_coefficient(c, bonds.in, -ω.b)
            end
        end
    end
    Random.seed!(1234)
    scenarios = [
        [(t, s[t]) for t in 1:5] for
        s in collect(Iterators.product([nothing], Ω, Ω, Ω, Ω))[:]
    ]
    probabilities = prod.(collect(Iterators.product([1.0], P, P, P, P))[:])
    SDDP.train(
        model;
        risk_measure = risk_measure,
        iteration_limit = 320,
        print_level = 0,
        sampling_scheme = CyclicHistorical(scenarios, 1),
    )
    Random.seed!(4321)
    simulations = SDDP.simulate(
        model,
        length(scenarios),
        [:stocks, :consumption];
        sampling_scheme = CyclicHistorical(scenarios, 1),
    )
    Z = [sum(stage[:stage_objective] for stage in s) for s in simulations]
    data = Dict{Float64,Float64}()
    for (z, p) in zip(Z, probabilities)
        if !haskey(data, z)
            data[z] = 0.0
        end
        data[z] += p
    end
    D = Distributions.DiscreteNonParametric(
        collect(keys(data)),
        collect(values(data)),
    )
    return (
        initial_stock = simulations[1][1][:stocks].out,
        consumption = Statistics.quantile(D, [0.0, 0.1, 0.5, 0.9, 1.0]),
    )
end

function build_figure_6()
    @info "Building Figure 6"
    γ_entropic = 0:1:50
    data_entropic = Dict(
        γ => d for
        (γ, d) in map(γ -> (γ, _run_instance(SDDP.Entropic(γ))), γ_entropic)
    )
    open(joinpath(@__DIR__, "figure_6_entropic.dat"), "w") do io
        return _write_matrix(
            io,
            hcat(
                γ_entropic,
                hcat([data_entropic[γ].consumption for γ in γ_entropic]...)',
            ),
        )
    end
    open(joinpath(@__DIR__, "figure_6_entropic_initial.dat"), "w") do io
        return _write_matrix(
            io,
            hcat(
                γ_entropic,
                reshape(
                    [data_entropic[γ].initial_stock for γ in γ_entropic],
                    (length(γ_entropic), 1),
                ),
            ),
        )
    end
    γ_avar = 0:0.01:1
    data_avar = Dict(
        γ => d for
        (γ, d) in map(γ -> (γ, _run_instance(SDDP.AVaR(1 - γ))), γ_avar)
    )
    open(joinpath(@__DIR__, "figure_6_avar.dat"), "w") do io
        return _write_matrix(
            io,
            hcat(γ_avar, hcat([data_avar[γ].consumption for γ in γ_avar]...)'),
        )
    end
    open(joinpath(@__DIR__, "figure_6_avar_initial.dat"), "w") do io
        return _write_matrix(
            io,
            hcat(
                γ_avar,
                reshape(
                    [data_avar[γ].initial_stock for γ in γ_avar],
                    (length(γ_avar), 1),
                ),
            ),
        )
    end
    return
end

build_figure_6()

###
### FIGURE 7
###

function _compute_stocks_bonds(risk_measure)
    Ω, P = [(s = 1.11, b = 1.02), (s = 1.04, b = 1.06)], [0.2, 0.8]
    scenarios = [
        [(t, s[t]) for t in 1:5] for
        s in collect(Iterators.product([(s = 1.0, b = 1.0)], Ω, Ω, Ω, Ω))[:]
    ]
    probabilities = prod.(collect(Iterators.product([1.0], P, P, P, P))[:])
    z_stocks_only = [prod(x[2].s for x in s) for s in scenarios]
    z_bonds_only = [prod(x[2].b for x in s) for s in scenarios]
    q = similar(probabilities)
    SDDP.adjust_probability(
        risk_measure,
        q,
        probabilities,
        SDDP.Noise{Any}[],
        z_stocks_only,
        false,
    )
    stocks_only = q' * z_stocks_only
    SDDP.adjust_probability(
        risk_measure,
        q,
        probabilities,
        SDDP.Noise{Any}[],
        z_bonds_only,
        false,
    )
    bonds_only = q' * z_bonds_only
    return stocks_only, bonds_only
end

function _compute_nested_avar(z, γ)
    function _eoh_cvar(z, p, γ)
        y = zeros(length(z))
        α = 0.0
        for k in sortperm(z; rev = false)
            if α >= (1 - γ)
                break
            end
            y[k] = min(p[k], (1 - γ) - α) / (1 - γ)
            α += y[k] * (1 - γ)
        end
        if sum(y) ≈ 0.0
            _, i = findmin(z)
            y[i] = 1.0
        end
        return z' * y
    end
    p = [0.2, 0.8]
    return _eoh_cvar(
        [
            _eoh_cvar(
                [
                    _eoh_cvar(
                        [
                            _eoh_cvar(z[1] * z[1] * z[1] .* z, p, γ)
                            _eoh_cvar(z[1] * z[1] * z[2] .* z, p, γ)
                        ],
                        p,
                        γ,
                    ),
                    _eoh_cvar(
                        [
                            _eoh_cvar(z[1] * z[2] * z[1] .* z, p, γ),
                            _eoh_cvar(z[1] * z[2] * z[2] .* z, p, γ),
                        ],
                        p,
                        γ,
                    ),
                ],
                p,
                γ,
            ),
            _eoh_cvar(
                [
                    _eoh_cvar(
                        [
                            _eoh_cvar(z[2] * z[1] * z[1] .* z, p, γ)
                            _eoh_cvar(z[2] * z[1] * z[2] .* z, p, γ)
                        ],
                        p,
                        γ,
                    ),
                    _eoh_cvar(
                        [
                            _eoh_cvar(z[2] * z[2] * z[1] .* z, p, γ),
                            _eoh_cvar(z[2] * z[2] * z[2] .* z, p, γ),
                        ],
                        p,
                        γ,
                    ),
                ],
                p,
                γ,
            ),
        ],
        p,
        γ,
    )
end

function build_figure_7()
    @info "Building Figure 7"
    open(joinpath(@__DIR__, "figure_7_entropic.dat"), "w") do io
        entropic = _compute_stocks_bonds.(SDDP.Entropic.(0:1:200))
        return _write_matrix(
            io,
            hcat(0:1:200, [e[1] for e in entropic], [e[2] for e in entropic]),
        )
    end
    open(joinpath(@__DIR__, "figure_7_avar.dat"), "w") do io
        eoh_avar = _compute_stocks_bonds.(SDDP.AVaR.(1.0:-0.01:0.0))
        return _write_matrix(
            io,
            hcat(
                0.0:0.01:1.0,
                [e[1] for e in eoh_avar],
                [e[2] for e in eoh_avar],
                _compute_nested_avar.(Ref([1.11, 1.04]), 0:0.01:1),
                _compute_nested_avar.(Ref([1.02, 1.06]), 0:0.01:1),
            ),
        )
    end
    return
end

build_figure_7()
