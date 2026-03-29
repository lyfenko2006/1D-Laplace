module Laplace_1D

using LinearAlgebra
export linear_solv, create_task

# решение СЛАУ методом прогонки
function linear_solv(A, b)
    if (size(A, 1) != size(A, 2))
        throw("Argument error")
    end
    n = size(A, 1)
    # прямой ход прогонки
    α = zeros(n)
    β = zeros(n)
    α[2] = -diag(A, 1)[1] / diag(A, 0)[1]
    β[2] = b[1] / diag(A, 0)[1]
    for i in 2:n - 1
        α[i + 1] = -diag(A, 1)[i] / (diag(A, -1)[i - 1] * α[i] + diag(A, 0)[i])
        β[i + 1] = (b[i] - diag(A, -1)[i - 1] * β[i]) / (diag(A, -1)[i - 1] * α[i] + diag(A, 0)[i])
    end
    # обратный ход прогонки
    y = zeros(n)
    y[n] = (b[n] - diag(A, -1)[n - 1] * β[n]) / (diag(A, -1)[n - 1] * α[n] + diag(A, 0)[n])
    for i in n - 1: -1 : 1
        y[i] = α[i + 1] * y[i + 1] + β[i + 1]
    end
    return y
end
# постановка задачи
function create_task(u, f, N)
    # Создали сетку
    x = [i / N for i in 0:N]
    # Создали правую часть СЛАУ
    b = zeros(N - 1)
    h = 1 / N
    b[1:N - 1] .+= -f.(x[2:N])
    b[1] += u(0) / h^2
    b[N - 1] += u(1) / h^2
    b .*= h^2
    # Создали матрицу СЛАУ
    dl = fill(-1, N - 2)
    du = fill(-1, N - 2)
    d = fill(2, N - 1)
    A = Tridiagonal(dl, d, du)
    return (A, b)
end

end # module Laplace_1D
