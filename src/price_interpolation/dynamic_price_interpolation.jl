#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

"""
    DynamicPriceInterpolation(; kwargs...)

Constuctor for the dynamic price interpolation value function described in
Downward, A., Dowson, O., and Baucke, R. (2018). On the convergence of a cutting
plane method for multistage stochastic programming problems with stagewise
dependent price uncertainty. Optimization Online.

### Keyword arguments
 - `dynamics`: a function that takes four arguments
        1. `price`: a Float64 (if uni-variate) or NTuple{N,Float64} if multi-variate
            that gives the price in the previous stage.
        2. `noise`: a single `NoiseRealization` of the price noise observed at
            the start of the stage.
        3. `t::Int`: the index of the stage of the problem t=1, 2, ..., T.
        4. `i::Int`: the markov state index of the problem i=1, 2, ..., S(t).
        The function should return a Float64 (if uni-variate) or NTuple{N,Float64}
        if multi-variate that gives the price for the current stage.
 - `initial_price`: a Float64 (if uni-variate) or NTuple{N,Float64} if multi-variate
    that gives the an initial value for each dimension of the price states.
 - `min_price`: a Float64 (if uni-variate) or NTuple{N,Float64} if multi-variate
    that gives the minimum value of each dimension of the price states.
 - `max_price`: a Float64 (if uni-variate) or NTuple{N,Float64} if multi-variate
    that gives the maximum value of each dimension of the price states.
 - `noise`: a finite-discrete distribution generated by `DiscreteDistribution`
 - `lipschitz_constant`: the maximum absolute subgradient of any price dimension
    in the domain bounded by `min_price` and `max_price`.

# Examples

A uni-variate price process:

    DynamicPriceInterpolation(
        dynamics = (price, noise, t, i) -> begin
                return price + noise - t
            end,
        initial_price = 50.0
        min_price = 0.0,
        max_price = 100.0,
        noise = DiscreteDistribution([-10.0, 40.0], [0.8, 0.2]),
        lipschitz_constant = 1e4
    )

A multi-variate price process:

    DynamicPriceInterpolation(
        dynamics = (price, noise, t, i) -> begin
                return (noise * price[1], noise * price[2] - noise)
            end,
        initial_price = (50.0, 50.0),
        min_price = (0.0,0.0),
        max_price = (100.0,100.0),
        noise = DiscreteDistribution([0.5, 1.2], [0.8, 0.2]),
        lipschitz_constant = 1e4
    )
"""
function DynamicPriceInterpolation(;
    dynamics = (p,w,t,i)->p,
    initial_price = 0.0,
    min_price = 0.0,
    max_price = 1.0,
    noise         = DiscreteDistribution([0.0]),
    lipschitz_constant = 1e6,
    cut_oracle = DefaultDynamicOracle(typeof(initial_price))
    )
    DynamicPriceInterpolation(
        initial_price,
        initial_price,
        min_price,
        max_price,
        noise,
        (p)->QuadExpr(p),
        dynamics,
        JuMP.Variable[],
        cut_oracle,
        lipschitz_constant,
        0.0
    )
end

# ==============================================================================

summarise{V<:DynamicPriceInterpolation}(::Type{V}) = "Dynamic Value Function"

# ==============================================================================

# stage, markov, price, cut
asynccutstoragetype{C,T,T2}(::Type{DynamicPriceInterpolation{C,T,T2}}) = Tuple{Int, Int, T, Cut}

# ==============================================================================

interpolate(vf::DynamicPriceInterpolation) = interpolate(vf.location, vf.mu)
function interpolate(price, mu)
    mu[1] + innerproduct(price, mu[2:end])
end
innerproduct(p::Float64, v::JuMP.Variable) = p * v
innerproduct(p::Float64, v::Vector{JuMP.Variable}) = (@assert length(v) == 1; p * v[1])
innerproduct{N,T}(p::NTuple{N,T}, v::Vector{JuMP.Variable}) = sum(p[i] * v[i] for i in 1:length(v))

# ==============================================================================

function _initializevaluefunction{V<:DynamicPriceInterpolation}(vf::V, m::JuMP.Model, sense, bound, N::Int)
    vf.bound = bound
    push!(vf.mu, @variable(m))
    for i in 1:N
        push!(vf.mu, @variable(m, lowerbound=-vf.lipschitz_constant, upperbound=vf.lipschitz_constant))
    end
    if 1 < N <= 4
        for price in Base.product(zip(vf.minprice, vf.maxprice)...)
            addpricecut!(sense, m, price, vf, bound)
        end
    else
        addpricecut!(sense, m, vf.minprice, vf, bound)
        addpricecut!(sense, m, vf.maxprice, vf, bound)
    end
    vf
end

# one-dimensional case
function initializevaluefunction{C,T<:Real, T2}(vf::DynamicPriceInterpolation{C,T, T2}, m::JuMP.Model, sense, bound)
    _initializevaluefunction(vf, m, sense, bound, 1)
end
# multi-dimensional case
function initializevaluefunction{C,N,T, T2}(vf::DynamicPriceInterpolation{C,NTuple{N,T}, T2}, m::JuMP.Model, sense, bound)
    _initializevaluefunction(vf, m, sense, bound, N)
end

# ==============================================================================

addpricecut!(::Max, sp, price, vf, affexpr) = @constraint(sp, interpolate(price, vf.mu) <= affexpr)
addpricecut!(::Min, sp, price, vf, affexpr) = @constraint(sp, interpolate(price, vf.mu) >= affexpr)

# ==============================================================================
#   updatevaluefunction!

function updatevaluefunction!{V<:DynamicPriceInterpolation}(m::SDDPModel{V}, settings::Settings, t::Int, sp::JuMP.Model)

    vf = valueoracle(sp)
    ex = ext(sp)

    current_price = vf.location

    # build the cut based on things stored in memory
    cut = constructcut(m, sp, ex, t, current_price)

    # if necessary, write to file
    if !settings.is_asyncronous && isopen(settings.cut_output_file)
        writecut!(settings.cut_output_file, ex.stage, ex.markovstate, current_price, cut)
    end

    # add the cut to oracle and JuMP model
    addcut!(m, sp, current_price, cut)

    # store this cut in m.ext[:cuts] for parallel if necessary
    if settings.is_asyncronous
        storeasynccut!(m, sp, current_price, cut)
    end
end

# ==============================================================================

function addcut!{V<:DynamicPriceInterpolation}(m::SDDPModel{V}, sp::JuMP.Model, current_price, cut::Cut)
    # get the value oracle
    vf = valueoracle(sp)

    # store the cut in the oracle
    storecut!(vf.oracle, m, sp, cut, current_price)

    # build affine expression
    affexpr = cuttoaffexpr(sp, cut)
    # add to JuMP model
    addpricecut!(ext(sp).sense, sp, current_price, vf, affexpr)

end

# add the cut from parallel
function addasynccut!{V<:DynamicPriceInterpolation,T}(m::SDDPModel{V}, cut::Tuple{Int, Int, T, Cut})
    sp = getsubproblem(m, cut[1], cut[2])
    addcut!(m, sp, cut[3], cut[4])
end

# ==============================================================================

rebuildsubproblem!(m::SDDPModel{V}, sp::JuMP.Model) where V <: DynamicPriceInterpolation = nothing

function rebuildsubproblem!(m::SDDPModel{DynamicPriceInterpolation{NanniciniOracle{T},T,T2}}, sp::JuMP.Model) where T where T2
    vf = valueoracle(sp)
    n = n_args(m.build!)
    ex = ext(sp)
    empty!(ex.states)
    empty!(ex.noises)

    sp2 = Model(solver = sp.solver)
    sp2.ext[:SDDP] = ex

    N = length(vf.mu) - 1
    empty!(vf.mu)
    _initializevaluefunction(vf, sp2, ext(sp2).sense, vf.bound, N)

    if n == 2
        m.build!(sp2, ex.stage)
    elseif n == 3
        m.build!(sp2, ex.stage, ex.markovstate)
    end
    for cut in validcuts(vf.oracle)
        affexpr = cuttoaffexpr(sp2, cut[1])
        addpricecut!(ext(sp2).sense, sp2, cut[2], vf, affexpr)
    end
    m.stages[ex.stage].subproblems[ex.markovstate] = sp2
end

# ==============================================================================

function processvaluefunctiondata{C,T2}(vf::DynamicPriceInterpolation{C,Float64, T2}, is_minimization::Bool, states::Union{Float64, AbstractVector{Float64}}...)
    cuts = [c[1] for c in vf.oracle.cuts]
    prices = [c[2] for c in vf.oracle.cuts]
    _processvaluefunctiondata(prices, cuts, vf.minprice, vf.maxprice, is_minimization, vf.lipschitz_constant, vf.bound, states...)
end

function postsolve!(::Type{ForwardPass}, m::SDDPModel{DynamicPriceInterpolation{V,T,T2}}, sp::JuMP.Model) where V <: NanniciniOracle where T where T2
    # we also need to overload the solve on the forward pass to record when
    # cuts are utilized. Alternatively, we could do this on the backward pass
    # to sample more widely.
    oracle = valueoracle(sp).oracle
    for i in 0:(oracle.cutsinmodel-1)
        if abs(sp.linconstrDuals[end-i]) > 1e-6
            oracle.iterations_since_last_active[end-i] = 0
        else
            oracle.iterations_since_last_active[end-i] += 1
        end
    end
end
