using Pkg
Pkg.activate(".")
try
    Pkg.update()
catch
    nothing
end
using CropGrowth
using Random, LinearAlgebra, StatsBase, Optimization, OptimizationBBO
using DataFrames, CSV
using Dates, ProgressMeter
using UnicodePlots, CairoMakie, ColorSchemes
