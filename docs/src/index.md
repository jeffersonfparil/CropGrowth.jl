# CropGrowth.jl

```@contents
Pages = ["index.md"]
Depth = 3
```

[CropGrowth.jl](https://github.com/jeffersonfparil/CropGrowth.jl) is a [Julia](https://julialang.org/) package for modelling crop growth curves using the [generalised logistic function](https://en.wikipedia.org/wiki/Generalised_logistic_function):

## Quickstart

```julia
using Pkg
Pkg.add(CropGrowth)
using CropGrowth
df = simulate()
df_out = fitgrowthmodels(df, verbose=true)
```

## Model

```math
y(t) = {A + {{K-A} \over {C + (Qe^{-Bt})^{1/v}}}}
```

where:

- ``y(t)``: biomass at time ``t`` (not part of the struct)
- ``A``: lower asymptote (initial or minimum biomass)
- ``K``: positively affects the upper asymptote. This is the final or maximum biomass if:
    + ``C = 1.00``, since:
    + ``y_{max} = A + (K-A)/C^{1/v}``, then 
    + ``y_{max} = A + K - A``, therefore: 
    + ``y_{max} = K``
- ``C``: negatively affects the final or maximum biomass
- ``Q``: negatively affects initial or minimum biomass
- ``e``: Euler's number (~2.71828)
- ``B``: growth rate
- ``v``: asymmetry parameter (``v ≥ 0``; small values: fast growth early; large values: fast growth later)

To solve for ``t`` at specific ``y``:

```math
t(y) = -{{1} \over {B}} \log {\left( { {{1} \over {Q}} \left( \left( {K - A} \over {y - A} \right)^v - C \right) } \right) }
```

## Installation

Please install the latest version of [Julia](https://julialang.org/) (we aim to update [CropGrowth.jl](https://github.com/jeffersonfparil/CropGrowth.jl) for the latest Julia release):

```shell
curl -fsSL https://install.julialang.org | sh
type -a julia
```

Install the [CropGrowth.jl](https://github.com/jeffersonfparil/CropGrowth.jl) library in Julia:

- Stable version:

```julia
using Pkg
Pkg.add(CropGrowth)
```

- Development version:

```julia
using Pkg
Pkg.add("https://github.com/jeffersonfparil/CropGrowth.jl")
```

Test the installation in Julia:

```julia
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

- The input `DataFrame` must contain the required columns specified in the global variable `REQUIRED_COLUMNS = ["entries", "sites", "replications", "growing_periods", "time_points"]`, as well as at least one additional trait column.
- If the `DataFrame` contains more than one trait column, only the first trait column will be used.
- Combinations with fewer than `min_t` time points will be skipped.
- The function uses a progress bar to indicate the fitting process if `verbose=true`.

## Examples

```julia
using CropGrowth, StatsBase
df = simulate(n_entries=5, n_sites=1, n_replications=1, n_growing_periods=1, n_time_points_per_growing_period=5, seed=123)
df_out_0, skipped_combinations_0 = fitgrowthmodels(df, show_plots=true)
df_out_C, skipped_combinations_C = fitgrowthmodels(df, C=Dict(:init=>1.0, :lower=>0.0, :upper=>100.0), show_plots=true)
df_out_Q, skipped_combinations_Q = fitgrowthmodels(df, Q=Dict(:init=>1.0, :lower=>0.0, :upper=>100.0), show_plots=true)
df_out_CQ, skipped_combinations_CQ = fitgrowthmodels(
    df,
    C=Dict(:init=>1.0, :lower=>0.0, :upper=>5.0),
    Q=Dict(:init=>1.0, :lower=>0.0, :upper=>5.0),
    show_plots=true
)
mean(df_out_0.R²)
mean(df_out_C.R²)
mean(df_out_Q.R²)
mean(df_out_CQ.R²)
```

## API

```@index
```

```@autodocs
Modules = [CropGrowth]
```
