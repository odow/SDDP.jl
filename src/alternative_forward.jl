#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""
    AlternativeForwardPass(
        forward_model::SDDP.PolicyGraph{T};
        forward_pass::AbstractForwardPass = DefaultForwardPass(),
    )

A forward pass that simulates using `forward_model`, which may be different to
the model used in the backwards pass.

When using this forward pass, you should almost always pass
[`SDDP.AlternativePostIterationCallback`](@ref) to the `post_iteration_callback`
argument of [`SDDP.train`](@ref).

This forward pass is most useful when the `forward_model` is non-convex and we
use a convex approximation of the model in the backward pass.

For example, in optimal power flow models, we can use an AC-OPF formulation as
the `forward_model` and a DC-OPF formulation as the backward model.

For more details see the paper:

Rosemberg, A., and Street, A., and Garcia, J.D., and Valladão, D.M., and Silva,
T., and Dowson, O. (2021). Assessing the cost of network simplifications in
long-term hydrothermal dispatch planning models. IEEE Transactions on
Sustainable Energy. 13(1), 196-206.
"""
struct AlternativeForwardPass{T} <: AbstractForwardPass
    model::PolicyGraph{T}
    forward_pass::AbstractForwardPass
    lock::ReentrantLock

    function AlternativeForwardPass(
        model::PolicyGraph{T};
        forward_pass::AbstractForwardPass = DefaultForwardPass(),
    ) where {T}
        return new{T}(model, forward_pass, ReentrantLock())
    end
end

function forward_pass(
    ::PolicyGraph{T},
    options::Options,
    pass::AlternativeForwardPass{T},
) where {T}
    # No need for locks here, delegate threadsafety to pass.forward_pass.
    return forward_pass(pass.model, options, pass.forward_pass)
end

"""
    AlternativePostIterationCallback(forward_model::PolicyGraph)

A post-iteration callback that should be used whenever
[`SDDP.AlternativeForwardPass`](@ref) is used.
"""
struct AlternativePostIterationCallback{T}
    model::PolicyGraph{T}
end

function (callback::AlternativePostIterationCallback)(result::IterationResult)
    # Only one thread is allowed to update the callback model at a time.
    callback.lock() do
        slave_update(callback.model, result)
        return
    end
    return
end
