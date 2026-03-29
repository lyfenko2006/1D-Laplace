module Func_Laplace_1D

using LinearAlgebra
export linear_solv, create_task

"""
    linear_solv(A, b) -> y

    Решение СЛАУ методом прогонки

    # Аргументы
    `A` - трехдиагональная матрица
    `b` - правая часть СЛАУ

    # Пример
    ```jldoctest
        julia> linear_solv(Tridiagonal([-1, 2], [3, 4, 3], [1, 2]), [1, -3, 5])
        3-element Vector{Float64}:
        1.0
        -2.0
        3.0000000000000004
    ```
"""
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

"""
    create_task(lower_bound, upper_bound, f, N) -> (A, b)

    Постановка обратной задачи одномерного уравнения Лапласа на отрезке [0, 1]

    # Аргументы
    `lower_bound` - нижнее граничное условие
    `upper_bound` - верхнее граничное условие
    `f` - вторая производная неизвестной функции
    `N` - число точек для разбиения отрезка [0, 1]

    # Пример
    ```jldoctest
        julia> N = 3
        julia> u = x -> sin(3x) * cos(4x)
        julia> f = x -> -25 * sin(3x) * cos(4x) - 24 * sin(4x) * cos(3x)
        julia> A, b = create_task(u(0), u(1), f, N)
        (Tridiagonal([-1], [2, 2], [-1]), [1.9502229726878755, -2.8459733453150107])
    ```
"""
function create_task(lower_bound, upper_bound, f, N)
    # Создали сетку
    x = [i / N for i in 0:N]
    # Создали правую часть СЛАУ
    b = zeros(N - 1)
    h = 1 / N
    b[1:N - 1] .+= -f.(x[2:N])
    b[1] += lower_bound / h^2
    b[N - 1] += upper_bound / h^2
    b .*= h^2
    # Создали матрицу СЛАУ
    dl = fill(-1, N - 2)
    du = fill(-1, N - 2)
    d = fill(2, N - 1)
    A = Tridiagonal(dl, d, du)
    return (A, b)
end

end # module Func_Laplace_1D