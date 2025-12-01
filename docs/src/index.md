# CropGrowth.jl

```@contents
Pages = ["index.md"]
Depth = 3
```

[CropGrowth.jl](https://github.com/jeffersonfparil/CropGrowth.jl) is a [Julia](https://julialang.org/) package for modelling crop growth curves using the [generalised logistic function](https://en.wikipedia.org/wiki/Generalised_logistic_function):

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
