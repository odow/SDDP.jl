#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
    return SDDP.train(
        convex;
        forward_pass = AlternativeForwardPass(nonconvex),
        parallel_scheme = AlternativeParallelScheme(nonconvex),
        kwargs...,
    )
end
