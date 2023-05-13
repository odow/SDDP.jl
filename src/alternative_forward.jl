#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""
    AlternativeForwardPass(;
        forward_model::SDDP.PolicyGraph{T},
    )

A forward pass for holding an alternative foward model. This is useful for
calculating cuts and simulating the policy under different models.
For example, simulating the forward passes with the AC-OPF to obtain forward trajectories,
and then refine the value functions on the backward pass using some convex approximation 
(e.g., DC with line losses) as in the following paper:
 - Rosemberg, Andrew Rosemberg, Alexandre Street, Joaquim Dias Garcia, Davi M. Valladão, Thuener Silva, and Oscar Dowson.
"Assessing the cost of network simplifications in long-term hydrothermal dispatch planning models." 
IEEE Transactions on Sustainable Energy 13, no. 1 (2021): 196-206.
"""
struct AlternativeForwardPass{T} <: AbstractForwardPass
    model::PolicyGraph{T}
end

function forward_pass(
    ::PolicyGraph{T},
    options::Options,
    pass::AlternativeForwardPass{T},
) where {T}
    return forward_pass(pass.model, options, DefaultForwardPass())
end

struct AlternativeParallelScheme{T} <: AbstractParallelScheme
    model::PolicyGraph{T}
end

Base.show(io::IO, ::AlternativeParallelScheme) = print(io, "alternative")

interrupt(::AlternativeParallelScheme) = nothing

function master_loop(
    scheme::AlternativeParallelScheme{T},
    model::PolicyGraph{T},
    options::Options,
) where {T}
    _initialize_solver(model; throw_error = false)
    while true
        result = iteration(model, options)
        slave_update(scheme.model, result)
        log_iteration(options)
        if result.has_converged
            return result.status
        end
    end
    return
end

function train_with_forward_model(nonconvex, convex; kwargs...)
    @assert isempty(setdiff(keys(nonconvex.initial_root_state), keys(convex.initial_root_state)))
    @assert isempty(setdiff(keys(convex.initial_root_state), keys(nonconvex.initial_root_state)))
    return SDDP.train(
        convex;
        forward_pass = AlternativeForwardPass(nonconvex),
        parallel_scheme = AlternativeParallelScheme(nonconvex),
        kwargs...,
    )
end
