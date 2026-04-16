module Func_validation

using Distances
export error_Euclidian, error_Chebyshev

"""
    error_Euclidian(aprox_val, foo, N) -> dist

    Расчет евклидова расстояния между значениями функции и ее апроксимации в N точках отрезка [0, 1]

    # Аргументы
    `aprox_val` - значения апроксимирующей функции на равномерной сетке
    `foo` - исходная функция
    `N` - количество точек для разбиения отрезка [0, 1]

    # Возвращаемое значение
    Евклидово расстояния между значениями функции и ее апроксимации в N точках отрезка [0, 1]

    # Пример
    ```jldoctest
        julia> error_Euclidian([1, 4 / 9], x -> x^2, 3)
        0.8888888888888888
        julia> error_Euclidian([1 / 9, 4 / 9], x -> x^2, 3)
        0.0
    ```
"""
function error_Euclidian(aprox_val, foo, N)
    x = [i / N for i in 0:N]
    val = foo.(x)[2 : N]
    Euclidean()(val, aprox_val)
end

"""
    error_Chebyshev(aprox_val, foo, N) -> dist

    Расчет чебышевского расстояния между значениями функции и ее апроксимации в N точках отрезка [0, 1]

    # Аргументы
    `aprox_val` - значения апроксимирующей функции на равномерной сетке
    `foo` - исходная функция
    `N` - количество точек для разбиения отрезка [0, 1]

    # Возвращаемое значение
    Евклидово расстояния между значениями функции и ее апроксимации в N точках отрезка [0, 1]

    # Пример
    ```jldoctest
        julia> error_Chebyshev([1, 4 / 9], x -> x^2, 3)
        0.8888888888888888
        julia> error_Chebyshev([1 / 9, 4 / 9], x -> x^2, 3)
        0.0
    ```
"""
function error_Chebyshev(aprox_val, foo, N)
    x = [i / N for i in 0:N]
    val = foo.(x)[2 : N]
    Chebyshev()(val, aprox_val)
end

end # module validation