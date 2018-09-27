struct Child
    # The child node.
    node::Node
    # The probability of transitioning to the child node from the parent.
    probability::Float64
end

struct Noise
    # The noise term.
    term  # TODO(odow): make this a concrete type?
    # The probability of sampling the noise term.
    probability::Float64
end

struct State
    # The incoming state variable in the subproblem.
    incoming::JuMP.VariableRef
    # The outgoing state variable in the subproblem.
    outgoing::JuMP.VariableRef
    # The "fishing" constraint so we can query the dual.
    fishing_dual  # TODO(odow): type this appropriately
end

struct Node
    # The JuMP subproblem.
    model::JuMP.Model
    # A vector of the child nodes.
    children::Vector{Child}
    # A vector of the discrete stagewise-independent noise terms.
    noise_terms::Vector{Noise}
    # A function parameterize(model::JuMP.Model, noise::T) that modifies the
    # JuMP model based on the observation of the noise.
    parameterize::Function  # TODO(odow): make this a concrete type?
    # A list of the state variables in the model.
    states::Vector{State}
    # A risk measure to apply over the children. If there are no children, this
    # should just default to expectation.
    risk_measure::AbstractRiskMeasure  # TODO(odow): make this a concrete type?
end

struct PolicyGraph
    # The value of the initial state variables.
    root_state::Vector{Float64}
    # Children of the root node.
    children::Vector{Child}
    # A risk measure to apply at the root node.
    risk_measure::AbstractRiskMeasure  # TODO(odow): make this a concrete type?
end
