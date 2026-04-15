module Laplace_1D
include("Func_Laplace_1D.jl")
include("Func_validation.jl")
include("Func_graphic.jl")
using .Func_Laplace_1D, .Func_validation, .Func_graphic
using Distances, Plots

N_arr = [5, 10, 20, 100, 1000, 10000]
euc_err = Vector{Float64}([])
cheb_err = Vector{Float64}([])
u = x -> sin(3x) * cos(4x)
f = x -> -25 * sin(3x) * cos(4x) - 24 * sin(4x) * cos(3x)
for N in N_arr
    A, b = create_task(u(0), u(1), f, N)
    y = linear_solv(A, b)
    push!(euc_err, error_Euclidian(y, u, N))
    push!(cheb_err, error_Chebyshev(y, u, N))
end
pattern = 1 ./ N_arr .^ 2
show_euc_error(N_arr, euc_err, pattern)
show_cheb_error(N_arr, cheb_err, pattern)
N = 20
A, b = create_task(u(0), u(1), f, N)
y = linear_solv(A, b)
x = [i / N for i in 0:N]
y_true = u.(x)
show_compare_approx(x, y, y_true)

end # module Laplace_1D