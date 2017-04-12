A scenario tree is a linked series of vectors

# Five stages in serial
ScenarioTree(
    stages = 5
)


# A markov chain of 5 stages and 2 markov states with uniform transition probabilities
ScenarioTree(
    stages        = 5,
    markov_states = 2
)

# A markov chain of 5 stages and 2 markov states with time invariant transition probabilities
ScenarioTree(
    stages        = 5,
    markov_states = 2,
    transition    = [0.4 0.6; 0.6 0.4]
)

# A markov chain of 5 stages and 2 markov states with time varying transition probabilities
ScenarioTree(
    stages        = 3,
    markov_states = 2,
    transition    = [ [0.4 0.6; 0.6 0.4], [0.3 0.7; 0.7 0.3] ]
)


struct Stages
    n::Int
end
struct MarkovStates
    n::Int
end

import Base: *
*(x::Int, ::Type{Stages}) = Stages(x)
*(x::Int, ::Type{MarkovStates}) = MarkovStates(x)
*(s::Stages, m::MarkovStates) = (s.x, m.x)
5Stages

5Stages * 2MarkovStates
