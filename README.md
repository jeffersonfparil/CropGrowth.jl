# CropGrowth

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jeffersonfparil.github.io/CropGrowth.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jeffersonfparil.github.io/CropGrowth.jl/dev/)
[![Build Status](https://github.com/jeffersonfparil/CropGrowth.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeffersonfparil/CropGrowth.jl/actions/workflows/CI.yml?query=branch%3Amain)

Modelling crop growth using logistic curves

## Quickstart

```julia
using Pkg
Pkg.add(CropGrowth)
using CropGrowth
df = simulate()
df_out = fitgrowthmodels(df, verbose=true)
```

## Model

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

## Examples

```julia
using CropGrowth, StatsBase
df = simulate(n_entries=5, n_sites=1, n_replications=1, n_growing_periods=1, n_time_points_per_growing_period=5, seed=123)
df_out_0 = fitgrowthmodels(df, show_plots=true)
df_out_C = fitgrowthmodels(df, C=Dict(:init=>1.0, :lower=>0.0, :upper=>100.0), show_plots=true)
df_out_Q = fitgrowthmodels(df, Q=Dict(:init=>1.0, :lower=>0.0, :upper=>100.0), show_plots=true)
df_out_CQ = fitgrowthmodels(
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