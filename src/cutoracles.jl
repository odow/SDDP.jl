abstract AbstractCutOracle

initialise(t::AbstractCutOracle) = error("""You must define an initialisation method for your cut oracle.""")

"""
    storecut!(oracle::AbstractCutOracle, cut::Cut)

    This function adds the cut to the CutOracle.
"""
storecut!(oracle::AbstractCutOracle, cut::Cut) = error("""
    You must define the function
        storecut!(oracle::$(typeof(oracle)), m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int, cut)
    that is overloaded for your oracle of type $(typeof(oracle)).
""")
# more expansive method that can also be overloaded
storecut!(oracle::AbstractCutOracle, m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int, cut) = storecut!(oracle, cut)

"""
    validcuts(oracle::AbstactCutOracle)
    This function returns an iterable list of all the valid cuts contained within the oracle.
"""
validcuts(oracle::AbstractCutOracle) = error("""You must define the function validcuts(oracle::$(typeof(oracle))) that is overloaded for your
    oracle of type $(typeof(oracle)).""")


struct DefaultCutOracle <: AbstractCutOracle
    cuts::Vector{Cut}
end
DefaultCutOracle() = DefaultCutOracle(Cut[])
initialise(::DefaultCutOracle, m, stage, markov, price) = DefaultCutOracle()
storecut!(oracle::DefaultCutOracle, m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int, cut) = push!(oracle.cuts, cut)
validcuts(oracle::DefaultCutOracle) = oracle.cuts
