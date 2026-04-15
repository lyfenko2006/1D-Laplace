module Func_graphic
using Plots

export show_euc_error, show_cheb_error, show_compare_approx

function show_euc_error(N_arr, euc_err, pattern)
    plot(N_arr, [euc_err pattern], lable = "Error")
    plot!(xscale=:log10, yscale=:log10)
    title!("Error in the Euclidean norm")
    xlabel!("N")
    ylabel!("Error")
    savefig("Euclidean_error.png")
end

function show_cheb_error(N_arr, cheb_err, pattern)
    plot(N_arr, [cheb_err pattern], lable = "Error")
    plot!(xscale=:log10, yscale=:log10)
    title!("Error in the Chebyshev norm")
    xlabel!("N")
    ylabel!("Error")
    savefig("Chebyshev_error.png")
end

function show_compare_approx(x, y, y_true)
    plot(x[2 : end - 1], [y y_true[2 : end - 1]], label = ["approx" "sin(3x)cos(4x)"])
    plot!(legend=:outerbottom, legendcolumns=2)
    title!("Compare functions")
    xlabel!("x")
    ylabel!("u")
    savefig("Compare_functions.png")
end

end #module Func_graphic