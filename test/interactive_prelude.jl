using Pkg
Pkg.activate(".")
try
    Pkg.update()
catch
    nothing
end
using Random, LinearAlgebra, StatsBase, Optimization, OptimizationBBO
using DataFrames, CSV
using Dates, ProgressMeter
using UnicodePlots, CairoMakie, ColorSchemes
