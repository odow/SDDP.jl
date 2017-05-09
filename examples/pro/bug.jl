module Foo

using JuMP, Gurobi

mutable struct T1
    x::Vector{Int}
end
T1() = T1(Int[])

function foo(f!::Function, n::Int)
    t1 = T1()
    for i in 1:n
        push!(t1.x, f!(i))
    end
    t1
end

end

using Foo

function serialize_round_trip(x)
    io = IOBuffer()
    serialize(io, x)
    seekstart(io)
    return deserialize(io)
end

@everywhere begin
    A = ones(10)
end

y = Foo.foo(3) do t
    # t2 = deepcopy(t)
    function bar(p)
        p * A[t]
    end
    A[t]
end

io = IOBuffer()
serialize(io, y)
seekstart(io)
@show take!(io)
@show serialize_round_trip(y)
