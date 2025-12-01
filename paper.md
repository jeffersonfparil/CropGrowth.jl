---
title: 'CropGrowth: Modelling crop growth using logistic curves in Julia'
tags:
  - Julia
  - crop growth modelling
  - agronomy
  - crop production
  - agriculture
authors:
  - name: Jefferson F. Paril
    orcid: 0000-0002-5693-4123
    affiliation: "1"
affiliations:
 - name: Agriculture Victoria, Bundoora, Victoria, Australia
   index: 1
date: 01 December 2025
bibliography: paper.bib
---

# Summary

`CropGrowth.jl` is a Julia package designed for simulating and fitting crop growth models using logistic curves. It provides tools for researchers and agronomists to analyze biomass accumulation over time, leveraging flexible parameterisation and statistical tools for model fitting.

# Statement of need

Understanding and predicting crop growth dynamics is critical for optimising agricultural productivity and ensuring food security. Traditional approaches to modelling crop growth often involve complex tools that may lack flexibility or require significant computational resources. `CropGrowth.jl` addresses these challenges by providing a lightweight, efficient, and user-friendly solution for modelling crop growth using logistic curves. 

This package is particularly valuable for agronomists, researchers, and students who need to simulate crop growth, fit models to experimental data, and visualize results. By leveraging Julia's high-performance capabilities [@Bezanson], `CropGrowth.jl` enables rapid prototyping and analysis, making it an essential tool for modern agricultural research.

# Mathematical model

$$
y(t) = {A + {{K-A} \over {C + (Qe^{-Bt})^{1/v}}}}
$$

where:

- $y(t)$: biomass at time $t$ (not part of the struct)
- $A$: lower asymptote (initial or minimum biomass)
- $K$: positively affects the upper asymptote. This be the final or maximum biomass if:
    + $C = 1.00$, since:
    + $y_{max} = A + (K-A)/C^{1/v}$, then 
    + $y_{max} = A + K - A$, therefore: 
    + $y_{max} = K$
- $C$: negatively affects the final or maximum biomass
- $Q$: negatively affects initial or minimum biomass
- $e$: Euler's number (~2.71828)
- $B$: growth rate
- $v$: asymmetry parameter ($v ≥ 0$; small values: fast growth early; large values: fast growth later)

To solve for $t$ at specific $y$:

$$
t(y) = -{{1} \over {B}} \log {\left( { {{1} \over {Q}} \left( \left( {K - A} \over {y - A} \right)^v - C \right) } \right) }
$$

# Usage and implementation details

## Quickstart

```julia
using Pkg
Pkg.add(CropGrowth)
using CropGrowth
df = simulate()
df_out = fitgrowthmodels(df, verbose=true)
```

## Implementation details

The main function, `fitgrowthmodels` fits generalized logistic growth models to the data provided in the input `DataFrame` and returns a `DataFrame` containing the fitted parameters, fit statistics, and time to reach specified percentages of the final value. The complete function signature is:

```julia
fitgrowthmodels(
  df::DataFrame;
  A = Dict(
      :init=>minimum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
      :lower=>0.0,
      :upper=>maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
  ),
  K = Dict(
      :init=>maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
      :lower=>0.0,
      :upper=>2*maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
  ),
  C = Dict(:init=>1.0, :lower=>1.0, :upper=>1.0),
  Q = Dict(:init=>1.0, :lower=>1.0, :upper=>1.0),
  B = Dict(:init=>1.0, :lower=>0.0, :upper=>10.0),
  v = Dict(:init=>1.0, :lower=>1e-5, :upper=>10.0),
  min_t::Int64 = 3,
  perc_of_final::Vector{Float64} = [0.5, 0.9],
  fit_statistic::String = "R²",
  maxiters::Int64 = 10_000,
  seed::Int64 = 42,
  show_plots::Bool = false,
  verbose::Bool = false,
)::Tuple{DataFrame, Vector{String}}
```

### Arguments

- `df::DataFrame`: Input data containing the required columns specified in `REQUIRED_COLUMNS = ["entries", "sites", "replications", "growing_periods", "time_points"]` and at least one trait column.
- `A::Dict`: Search space for the parameter `A` (lower asymptote). Contains `:init`, `:lower`, and `:upper` keys. Defaults to the minimum and maximum of the trait column with `init=minimum`.
- `K::Dict`: Search space for the parameter `K` (upper asymptote). Contains `:init`, `:lower`, and `:upper` keys. Defaults to the minimum and 2×maximum of the trait column with `init=maximum`.
- `C::Dict`: Search space for the parameter `C`. Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=1.0`, `upper=1.0`.
- `Q::Dict`: Search space for the parameter `Q`. Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=1.0`, `upper=1.0`.
- `B::Dict`: Search space for the parameter `B` (growth rate). Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=0.0`, `upper=10.0`.
- `v::Dict`: Search space for the parameter `v` (asymmetry parameter). Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=1e-5`, `upper=10.0`.
- `min_t::Int64`: Minimum number of time points required to fit the growth model for a specific combination of entry, site, replication, and growing period. Defaults to `3`.
- `perc_of_final::Vector{Float64}`: Percentages of the final value for which the time to reach these percentages will be calculated. Defaults to `[0.5, 0.9]`.
- `fit_statistic::String`: The fit statistic to be used for evaluating the model. Must be one of `["R²", "RMSE", "MSE", "MAE", "ρ"]`. Defaults to `"R²"`.
- `maxiters::Int64`: Maximum number of iterations allowed for the optimization process. Defaults to `10_000`.
- `seed::Int64`: Random seed for reproducibility. Defaults to `42`.
- `show_plots::Bool`: Whether to show fitted growth curve plots. Defaults to `false`.
- `verbose::Bool`: Whether to display progress and additional information during the fitting process. Defaults to `false`.

### Returns

`Tuple{DataFrame, Vector{String}}`: 
  + The first element is a `DataFrame` containing the fitted parameters (`A`, `K`, `C`, `Q`, `B`, `v`), fit statistics, value of the growth models at ``t=0`` (`y_t0`), maximum value of the growth model (`y_max`), and time to reach specified percentages of the final value for each combination of entry, site, replication, and growing period.
  + The second element is a `Vector{String}` containing the combinations that were skipped due to insufficient data points.

### Notes

- This package leverages `Optimization.jl` [@Vaibhav] to fit the logistic models, `DataFrames.jl` [@Bouchet-Valat] for efficient data handling, `CSV.jl` [@JuliaData] for file input and output, `ProgressMeter.jl` [@Holy] for status monitoring, and `UnicodePlots.jl` [@Danisch] for simple plotting.
- The input `DataFrame` must contain the required columns specified in the global variable `REQUIRED_COLUMNS`, as well as at least one additional trait column.
- If the `DataFrame` contains more than one trait column, only the first trait column will be used.
- Combinations with fewer than `min_t` time points will be skipped, and a warning will be issued.
- The function uses a progress bar to indicate the fitting process if `verbose=true`.


# Conclusion

`CropGrowth.jl` simplifies the modelling of crop growth dynamics, making it a valuable tool for agricultural research. For more details, visit the [documentation](https://jeffersonfparil.github.io/CropGrowth.jl/stable/).

# Acknowledgements

We acknowledge the developers and contributors of the packages used in `CropGrowth.jl`.
The complete list is found in the Project.toml of the repository.

# References