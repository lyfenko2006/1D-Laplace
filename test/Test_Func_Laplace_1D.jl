module Test_Func_Laplace_1D
include(Func_Laplace_1D.jl)
using .Func_Laplace_1D
using Test
@test linear_solv(Tridiagonal([-1, 2], [3, 4, 3], [1, 2]), [1, -3, 5]) ≈ [1, -2, 3]


end # module Test_Func_Laplace_1D