# module Test
#
#     struct T2
#         x::Int
#     end
#
#     struct T1
#         x::Vector{T2}
#     end
#
#     function T1(f::Function, n::Int)
#         y = T2[]
#         for i in 1:n
#             push!(y, f(i))
#         end
#         T1(y)
#     end

# end
#
# a = [1,2,3]
#
# y = Test.T1(3) do i
#     y = Test.T2(2)
#     function foo(x)
#         x*a[i] + y.x
#     end
#     y
# end
#
# @show Test.serialize_round_trip(y)

function serialize_round_trip(x)
    io = IOBuffer()
    serialize(io, x)
    seekstart(io)
    return deserialize(io)
end
global x = ones(10)
A = B = x

serialize_round_trip(A)

# @show Base.Serializer.handle_deserialize(ss, i)
# @show @which deserialize(ss)
