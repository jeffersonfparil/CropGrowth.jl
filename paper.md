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
affiliations:
 - name: Agriculture Victoria, AgriBio, Centre for AgriBioscience, 5 Ring Rd, Bundoora, VIC, Australia
   index: 1
date: 01 December 2025
bibliography: paper.bib
---

# Summary

`CropGrowth.jl` is a Julia package for modelling crop growth using the generalised logistic function. The package offers an easy-to-use interface, flexible parameterisation, and when combined with Julia’s speed and interactive environment, it provides an efficient and accessible tool for researchers, agronomists, and students working on crop growth analysis.

# Statement of need

Understanding and predicting crop growth dynamics is critical for optimising agricultural productivity and ensuring food security. Conventional approaches to crop growth modelling typically rely on numerous disparate tools, compelling users to assemble *ad hoc* solutions that are neither flexible nor intuitive. `CropGrowth.jl` addresses these challenges by providing a lightweight, efficient, and user-friendly solution for modelling crop growth using the generalised logistic function [@Richards]. It supports flexible parameterisation, from the full six-parameter model to simpler forms such as Gompertz function [@Gompertz] and the standard logistic curve.

This package offers significant utility for agronomists, researchers, and students seeking to characterise crop growth responses across one or more genotypes using field trial data, which may encompass multiple years, environments, genotypes, and treatment factors. By leveraging Julia's high-performance capabilities [@Bezanson], `CropGrowth.jl` enables rapid prototyping and analysis, making it an essential tool for modern agricultural research.

This package leverages `OptimizationBBO.jl` [@Vaibhav] to fit the logistic models (specifically the `BBO_adaptive_de_rand_1_bin_radiuslimited()` optimiser), `DataFrames.jl` [@Bouchet-Valat] for efficient data handling, `CSV.jl` [@Quinn] for table reading and writing, `ProgressMeter.jl` [@Holy] for status monitoring, and `UnicodePlots.jl` [@Stocker] for simple plotting.

# Mathematical model

$$
y(t) = {A + {{K-A} \over {C + (Qe^{-Bt})^{1/v}}}}
$$

where:

- $y(t)$: biomass at time $t$
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

`CropGrowth.jl` simplifies the modelling of crop growth dynamics, making it a valuable tool for agricultural research. For more details, visit the documentations ([stable version](https://jeffersonfparil.github.io/CropGrowth.jl/stable/) or [development version](https://jeffersonfparil.github.io/CropGrowth.jl/dev/)).

# Acknowledgements

We acknowledge Agriculture Victoria Research, Dairy Australia, and the Gardiner Foundation for their financial support, which made the development of this software possible. We also extend our appreciation to the developers and contributors of the packages used in `CropGrowth.jl`. The complete list of these dependencies is found in the [Project.toml](https://github.com/jeffersonfparil/CropGrowth.jl/blob/main/Project.toml) file within the repository.

# References