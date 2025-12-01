---
title: 'CropGrowth: Modelling crop growth using logistic curves in Julia'
tags:
  - Julia
  - crop growth modelling
  - agronomy
  - agriculture
authors:
  - name: Jefferson F. Paril
    orcid: 0000-0002-5693-4123
    affiliation: "1"
  - name: Noel O. I. Cogan
    orcid: 0000-0003-1074-8998
    affiliation: "1,2"
  - name: M. Michelle Malmberg
    orcid: 0000-0002-6351-5904
    affiliation: "1"
affiliations:
 - name: Agriculture Victoria, AgriBio, Centre for AgriBioscience, 5 Ring Rd, Bundoora, VIC, Australia
   index: 1
 - name: School of Applied Systems Biology, LaTrobe University, Bundoora, VIC, Australia
   index: 2
date: 01 December 2025
bibliography: paper.bib
---

# Summary

`CropGrowth.jl` is a Julia package designed for simulating and fitting crop growth models using logistic curves. It provides tools for researchers and agronomists to analyze biomass accumulation over time, leveraging flexible parameterisation and statistical tools for model fitting.

# Statement of need

Understanding and predicting crop growth dynamics is critical for optimising agricultural productivity and ensuring food security. Conventional approaches to crop growth modelling typically rely on numerous disparate tools, compelling users to assemble *ad hoc* solutions that are neither flexible nor intuitive. `CropGrowth.jl` addresses these challenges by providing a lightweight, efficient, and user-friendly solution for modelling crop growth using logistic curves. 

This package offers significant utility for agronomists, researchers, and students seeking to characterise crop growth responses across one or more genotypes using field trial data, which may encompass multiple years, environments, genotypes, and treatment factors. By leveraging Julia's high-performance capabilities [@Bezanson], `CropGrowth.jl` enables rapid prototyping and analysis, making it an essential tool for modern agricultural research.

This package leverages `Optimization.jl` [@Vaibhav] to fit the logistic models, `DataFrames.jl` [@Bouchet-Valat] for efficient data handling, `CSV.jl` [@JuliaDataContributors] for file input and output, `ProgressMeter.jl` [@HolyTim] for status monitoring, and `UnicodePlots.jl` [@PlotsContributors] for simple plotting.

# Mathematical model

$$
y(t) = {A + {{K-A} \over {C + (Qe^{-Bt})^{1/v}}}}
$$

where:

- $y(t)$: biomass at time $t$ (not part of the struct)
- $A$: lower asymptote (initial or minimum biomass)
- $K$: positively affects the upper asymptote. This is the upper asymptote or maximum biomass if $C = 1.00$.
- $C$: negatively affects the final or maximum biomass
- $Q$: negatively affects initial or minimum biomass
- $e$: Euler's number (~2.71828)
- $B$: growth rate
- $v$: asymmetry parameter ($v ≥ 0$; small values: fast growth early; large values: fast growth later)

To solve for $t$ at specific $y$:

$$
t(y) = -{{1} \over {B}} \log {\left( { {{1} \over {Q}} \left( \left( {K - A} \over {y - A} \right)^v - C \right) } \right) }
$$

# Conclusion

`CropGrowth.jl` simplifies the modelling of crop growth dynamics, making it a valuable tool for agricultural research. For more details, visit the [documentation](https://jeffersonfparil.github.io/CropGrowth.jl/stable/).

# Acknowledgements

We acknowledge the developers and contributors of the packages used in `CropGrowth.jl`. The complete list is found in the [Project.toml](https://github.com/jeffersonfparil/CropGrowth.jl/blob/main/Project.toml) of the repository.

# References