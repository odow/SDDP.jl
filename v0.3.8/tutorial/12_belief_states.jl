# # Advanced II: belief states

# `SDDP.jl` includes an implementation of the algorithm described in Dowson, O.,
# Morton, D.P., & Pagnoncelli, B. (2019). Partially observable multistage
# stochastic programming. [[link]](http://www.optimization-online.org/DB_HTML/2019/03/7141.html)

# Proper documentation will be forthcoming.

# In the mean-time, here is most of what you need to know.

# ## Defining the ambiguity partition

# Given a [`SDDP.Graph`](@ref) object (see [Create a general policy graph](@ref)
# for details), we can define the ambiguity partition using
# [`SDDP.add_ambiguity_set`](@ref).

# For example, first we create a Markovian graph:

using SDDP

G = SDDP.MarkovianGraph([[0.5 0.5], [0.2 0.8; 0.8 0.2]])

# Then we add  an ambiguity set over the nodes in the first stage:

SDDP.add_ambiguity_set(G, [(1, 1), (1, 2)])

# Then we add  an ambiguity set over the nodes in the second stage:

SDDP.add_ambiguity_set(G, [(2, 1), (2, 2)])

# This results in the graph:

G
