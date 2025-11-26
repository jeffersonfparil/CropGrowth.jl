# CropGrowth.jl

CropGrowth.jl is a Julia package for modelling crop growth curves using the generalised logistic function:

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
t(y) = -{{1} \over {B}} \log {\left( { \left( {K - A} \over {y - A} \right)^v - C } \over {Q} \right) }
$$


The `GenomicBreeding` module provides a comprehensive suite of tools for genomic prediction, genome-wide association studies (GWAS), and data handling in genomic breeding. It integrates functionalities from `GenomicBreedingCore`, `GenomicBreedingIO`, `GenomicBreedingModels`, and `GenomicBreedingPlots` to offer efficient and scalable solutions for genetic data analysis and visualisation.

```@contents
Pages = ["index.md"]
Depth = 3
```

## Installation

We designed [GenomicBreeding.jl](https://github.com/GenomicBreeding/GenomicBreeding.jl) to work on an HPC running Linux (the various components, i.e. [GenomicBreedingCore.jl](https://github.com/GenomicBreeding/GenomicBreedingCore.jl), [GenomicBreedingIO.jl](https://github.com/GenomicBreeding/GenomicBreedingIO.jl), [GenomicBreedingModels.jl](https://github.com/GenomicBreeding/GenomicBreedingModels.jl), and [GenomicBreedingPlots.jl](https://github.com/GenomicBreeding/GenomicBreedingPlots.jl) work on a single Linux PC too).

First, if you have not yet, please install [Julia](https://julialang.org/) (the latest version as we aim to update this package for the latest Julia release):

```shell
curl -fsSL https://install.julialang.org | sh
type -a julia
```

Install the [CropGrowth.jl](https://github.com/jeffersonfparil/CropGrowth.jl) library in Julia:

```julia
using Pkg
Pkg.add("https://github.com/jeffersonfparil/CropGrowth.jl")
```

## Examples

The main API is the `fitgrowthmodels()` function:

```julia
using Pkg
Pkg.add(CropGrowth)
using CropGrowth
df = simulate()
df_out = fitgrowthmodels(df, verbose=true)
```

## API

```@index
```

```@autodocs
Modules = [CropGrowth]
```
