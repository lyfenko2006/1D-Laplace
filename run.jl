import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include("src/Laplace_1D.jl")