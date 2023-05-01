# The function `_lattice_approximation` is derived from a function of the same name in the
# `ScenTrees.jl` package by Kipngeno Kirui and released under the MIT license.
# The reproduced function, and other functions contained only in this file, are also
# released under MIT.
#
# Copyright (c) 2019 Kipngeno Kirui <kipngenokirui1993@gmail.com>
# Copyright (c) 2019 Oscar Dowson <o.dowson@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

function find_min(x::Vector{T}, y::T) where {T<:Real}
    best_i = 0
    best_z = Inf
    for i in 1:length(x)
        z = abs(x[i] - y)
        if z < best_z
            best_i = i
            best_z = z
        end
    end
    return best_z, best_i
end

function _lattice_approximation(
    f::Function,
    states::Vector{Int},
    scenarios::Int,
)
    return _lattice_approximation(
        f,
        states,
        scenarios,
        [f()::Vector{Float64} for _ in 1:scenarios],
    )
end

function _quantiles(x, N)
    if N == 1
        return [Statistics.mean(x)]
    end
    return Statistics.quantile(x, range(0.01, 0.99; length = N))
end

function _lattice_approximation(
    f::Function,
    states::Vector{Int},
    scenarios::Int,
    simulations::Vector{Vector{Float64}},
)
    simulation_matrix = reduce(hcat, simulations)
    support = map(1:length(states)) do t
        return _quantiles(@view(simulation_matrix[t, :]), states[t])
    end
    probability = [zeros(states[t-1], states[t]) for t in 2:length(states)]
    prepend!(probability, Ref(zeros(1, states[1])))
    distance = 0.0
    for (n, path) in enumerate(simulations)
        dist, last_index = 0.0, 1
        for t in 1:length(states)
            for i in 1:length(states[t])
                if sum(@view probability[t][:, i]) < 1.3 * sqrt(n) / states[t]
                    support[t][i] = path[t]
                end
            end
            min_dist, best_idx = find_min(support[t], path[t])
            dist += min_dist^2
            probability[t][last_index, best_idx] += 1.0
            support[t][best_idx] -=
                min_dist * (support[t][best_idx] - path[t]) / (3000 + n)^0.75
            last_index = best_idx
        end
        distance = (distance * (n - 1) + dist) / n
    end
    for p in probability
        p ./= sum(p, dims = 2)
        if any(isnan, p)
            p[vec(isnan.(sum(p, dims = 2))), :] .= 0.0
        end
    end
    return support, probability
end

"""
    _allocate_support_budget(f, budget, scenarios)

Allocate the `budget` nodes amongst the stages for a Markovian approximation.
By default, we distribute nodes based on the relative variance of the stages.
"""
function _allocate_support_budget(
    f::Function,
    budget::Int,
    scenarios::Int,
)::Vector{Int}
    return _allocate_support_budget(
        [f()::Vector{Float64} for _ in 1:scenarios],
        budget,
        scenarios,
    )
end

function _allocate_support_budget(
    simulations::Vector{Vector{Float64}},
    budget::Int,
    scenarios::Int,
)::Vector{Int}
    stage_var = Statistics.var(simulations)
    states = ones(Int, length(stage_var))
    if budget < length(stage_var)
        @warn(
            "Budget for nodes is less than the number of stages. Using one " *
            "node per stage.",
        )
        return states
    end
    s = sum(stage_var)
    if s ≈ 0.0
        # If the sum of the variances is 0, then the simulator must be
        # deterministic. Regardless of the budget, return a single Markov state
        # for each stage.
        return states
    end
    for i in 1:length(states)
        states[i] = max(1, round(Int, budget * stage_var[i] / s))
    end
    while sum(states) != budget
        if sum(states) > budget
            states[argmax(states)] -= 1
        else
            states[argmin(states)] += 1
        end
    end
    return states
end

_allocate_support_budget(::Any, budget::Vector{Int}, ::Int) = budget

"""
    MarkovianGraph(
        simulator::Function;
        budget::Union{Int,Vector{Int}},
        scenarios::Int = 1000,
    )

Construct a Markovian graph by fitting Markov chain to scenarios generated by
`simulator()`.

`budget` is the total number of nodes in the resulting Markov chain. This can
either be specified as a single `Int`, in which case we will attempt to
intelligently distributed the nodes between stages. Alternatively, `budget` can
be a `Vector{Int}`, which details the number of Markov state to have in each
stage.
"""
function MarkovianGraph(
    simulator::Function;
    budget::Union{Int,Vector{Int}},
    scenarios::Int = 1000,
)
    scenarios = max(scenarios, 10)
    simulations = [simulator()::Vector{Float64} for _ in 1:scenarios]
    states = _allocate_support_budget(simulations, budget, scenarios)
    support, probability =
        _lattice_approximation(simulator, states, scenarios, simulations)
    g = Graph((0, 0.0))
    for (i, si) in enumerate(support[1])
        add_node(g, (1, si))
        add_edge(g, (0, 0.0) => (1, si), probability[1][1, i])
    end
    for t in 2:length(support)
        for (j, sj) in enumerate(support[t])
            add_node(g, (t, sj))
            for (i, si) in enumerate(support[t-1])
                add_edge(g, (t - 1, si) => (t, sj), probability[t][i, j])
            end
        end
    end
    return g
end
