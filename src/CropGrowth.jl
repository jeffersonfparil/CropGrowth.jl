module CropGrowth

using Random, LinearAlgebra, StatsBase, Optimization, OptimizationBBO
using DataFrames, CSV
using Dates, ProgressMeter
using UnicodePlots, CairoMakie, ColorSchemes

const REQUIRED_COLUMNS =
    ["entries", "sites", "replications", "growing_periods", "time_points"]
include("io.jl")
include("model.jl")
include("fit.jl")

export GrowthModel, fitstatistics, generalisedlogistic, modelgrowth, timetomaxperc
export REQUIRED_COLUMNS, simulate, readdelimited, writedelimited
export fitgrowthmodels

end
