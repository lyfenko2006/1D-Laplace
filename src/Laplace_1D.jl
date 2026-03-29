module Laplace_1D
include("Func_Laplace_1D.jl")
using .Func_Laplace_1D

N = 5
u = x -> sin(3x) * cos(4x)
f = x -> -25 * sin(3x) * cos(4x) - 24 * sin(4x) * cos(3x)
A, b = create_task(u(0), u(1), f, N)
y = linear_solv(A, b)

end # module Laplace_1D
