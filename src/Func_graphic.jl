module Func_graphic
using Plots

export show_euc_error, show_cheb_error, show_compare_approx

"""
    show_euc_error(N_arr, euc_err, pattern, path)

    Создание графика ошибки в Евклидовой норме и сравнение его с эталонным графиком
    
    # Аргументы
    N_arr - массив точек отрезка, соответствующих значениям функций
    euc_err - массив значений ошибки в данных точках
    pattern - массив с эталонными значениями ошибки (с точностью до константы) для данного типа сходимости
    path - путь для сохранения графика

"""

function show_euc_error(N_arr, euc_err, pattern, path)
    plot(N_arr, [euc_err pattern], label = ["Euclidean error" "1 / N^2"])
    plot!(xscale=:log10, yscale=:log10)
    title!("Error in the Euclidean norm")
    xlabel!("N")
    ylabel!("Error")
    savefig(joinpath(path, "Euclidean_error.png"))
end

"""
    show_cheb_error(N_arr, euc_err, pattern, path)

    Создание графика ошибки в Чебышевской норме и сравнение его с эталонным графиком
    
    # Аргументы
    N_arr - массив точек отрезка, соответствующих значениям функций
    cheb_err - массив значений ошибки в данных точках
    pattern - массив с эталонными значениями ошибки (с точностью до константы) для данного типа сходимости
    path - путь для сохранения графика
    
"""

function show_cheb_error(N_arr, cheb_err, pattern, path)
    plot(N_arr, [cheb_err pattern], label = ["Chebyshev_error" "1 / N^2"])
    plot!(xscale=:log10, yscale=:log10)
    title!("Error in the Chebyshev norm")
    xlabel!("N")
    ylabel!("Error")
    savefig(joinpath(path, "Chebyshev_error.png"))
end

"""
    show_compare_approx(x, y, y_true, path)

    Создание графика-сравнения точности апроксимации исходной функции данной
    
    # Аргументы
    x - массив точек отрезка, соответствующих значениям функций
    y - массив значений апроксимирующей функции в данных точках
    y_true - массив значений исходной функции в данных точках
    path - путь для сохранения графика
    
"""

function show_compare_approx(x, y, y_true, path)
    plot(x, [y y_true], label = ["approx" "sin(3x)cos(4x)"])
    # plot!(legend=:outerbottom, legendcolumns=2)
    title!("Compare functions")
    xlabel!("x")
    ylabel!("u")
    savefig(joinpath(path, "Compare_functions.png"))
end

end #module Func_graphic