"""
    struct GrowthModel

Holds parameters for the generalized logistic growth model:

```math
y(t) = A + \\frac{K-A}{C + (Qe^{-Bt})^{1/v}}
```

where:

- ``y(t)``: biomass at time ``t`` (not part of the struct)
- ``A``: lower asymptote (initial or minimum biomass)
- ``K``: positively affects the upper asymptote (can be the final or maximum biomass if ``C = 1.00``, since ``y_{max} = A + (K-A)/C^(1/v)``, then ``y_{max} = A + K - A``, therefore: ``y_{max} = K``)
- ``C``: negatively affects the final or maximum biomass
- ``Q``: negatively affects initial or minimum biomass
- ``e``: Euler's number (~2.71828)
- ``B``: growth rate
- ``v``: asymmetry parameter (``v ≥ 0``; small values: fast growth early; large values: fast growth later)

Additional information are provided in the struct:

- `y_max`: maximum value of the growth model (``y_{max} = A + (K-A)/C^(1/v)``)
- `fit_statistics`: a dictionary containing fit statistics such as R², RMSE, MSE, MAE, and Pearson's correlation coefficient (ρ)

# Constructor

```julia
GrowthModel(;
    A::Float64,
    K::Float64,
    C::Float64,
    Q::Float64,
    B::Float64,
    v::Float64,
    fit_statistics::Dict{String,Float64} = Dict(""=>NaN),
    ϵ::Float64 = 1e-12,
)
```

# Notes:
- The struct includes a constructor that automatically calculates `y_max` based on the provided parameters.
- The `fit_statistics` field is optional and defaults to a dictionary with a single entry of an empty string key and `NaN` value.
- A small epsilon value (`ϵ`) is added to the calculation of `y_max` to prevent division by zero errors.

# Example
```jldoctest; setup = :(using CropGrowth)
julia> growth_model = GrowthModel(A=0.0, K=10.0, C=1.0, Q=1.0, B=0.75, v=0.1);

julia> (growth_model.y_max - growth_model.K) < 1e-9
true

julia> (growth_model.y_max - 10.0) < 1e-9
true
```
"""
mutable struct GrowthModel
    A::Float64
    K::Float64
    C::Float64
    Q::Float64
    B::Float64
    v::Float64
    y_max::Float64
    fit_statistics::Dict{String,Float64}
    function GrowthModel(;
        A::Float64,
        K::Float64,
        C::Float64,
        Q::Float64,
        B::Float64,
        v::Float64,
        fit_statistics::Dict{String,Float64} = Dict(""=>NaN),
        ϵ::Float64 = 1e-12,
    )
        y_max = A + (K-A)/(C^(1/v) + ϵ)
        new(A, K, C, Q, B, v, y_max, fit_statistics)
    end
end

"""
    fitstatistics(; y::Vector{Float64}, ŷ::Vector{Float64})::Dict{String, Float64}

Compute various fit statistics for a given set of observed and predicted values.

# Arguments
- `y::Vector{Float64}`: Observed data points.
- `ŷ::Vector{Float64}`: Predicted data points.

# Returns
A dictionary containing the following fit statistics:
- `"R²"`: Coefficient of determination, a measure of how well the predicted values explain the variance in the observed data.
- `"rmse"`: Root mean squared error, a measure of the average magnitude of the residuals.
- `"mse"`: Mean squared error, the average of the squared residuals.
- `"mae"`: Mean absolute error, the average of the absolute residuals.
- `"ρ"`: Pearson correlation coefficient between the observed and predicted values.

# Notes
- Ensure that the `y` and `ŷ` vectors are of the same length and correspond to the same observations.
- The function assumes that both `y` and `ŷ` are non-empty vectors.

# Example
```jldoctest; setup = :(using CropGrowth)
julia> y1::Vector{Float64} = rand(10); y2::Vector{Float64} = rand(10);

julia> stats0 = fitstatistics(y=y1, ŷ=y1);

julia> stats1 = fitstatistics(y=y1, ŷ=y2);

julia> stats0["R²"] > stats1["R²"]
true

julia> stats0["rmse"] < stats1["rmse"]
true

julia> stats0["mse"] < stats1["mse"]
true

julia> stats0["mae"] < stats1["mae"]
true

julia> stats0["ρ"] > stats1["ρ"]
true
```
"""
function fitstatistics(; y::Vector{Float64}, ŷ::Vector{Float64})::Dict{String,Float64}
    # y::Vector{Float64} = rand(10); ŷ::Vector{Float64} = rand(10)
    if length(y) != length(ŷ)
        throw(ArgumentError("Vectors y and ŷ must be of the same length."))
    end
    ss_total = sum((y .- mean(y)) .^ 2)
    ss_residual = sum((y .- ŷ) .^ 2)
    R² = 1.0 - (ss_residual / ss_total)
    mae = mean(abs.(y .- ŷ))
    mse = mean((y .- ŷ) .^ 2)
    rmse = sqrt(mse)
    ρ = cor(y, ŷ)
    return Dict("R²" => R², "rmse" => rmse, "mse" => mse, "mae" => mae, "ρ" => ρ)
end

"""
    generalisedlogistic(t::Vector{Float64}; A=0.0, K=1.0, C=1.0, Q=1.0, B=1.0, v=1.0)::Vector{Float64}

Computes the generalized logistic function for a given vector of time points `t`.

# Arguments
- `t::Vector{Float64}`: A vector of time points at which to evaluate the function.

# Keyword Arguments
- `A::Float64` (default: `0.0`): lower asymptote (initial or minimum biomass)
- `K::Float64` (default: `1.0`): upper asymptote (can be the final or maximum biomass if ``C = 1.00``, since ``y_{max} = A + (K-A)/C^(1/v)``, then ``y_{max} = A + K - A``, therefore: ``y_{max} = K``)
- `C::Float64` (default: `1.0`): negatively affects the final or maximum biomass
- `Q::Float64` (default: `1.0`): negatively affects initial or minimum biomass
- `B::Float64` (default: `1.0`): growth rate
- `v::Float64` (default: `1.0`): asymmetry parameter (``v ≥ 0``; small values: fast growth early; large values: fast growth later)

# Returns
- `Vector{Float64}`: A vector of values representing the generalized logistic function evaluated at each time point in `t`.

# Example
```jldoctest; setup = :(using CropGrowth)
julia> t = collect(0.0:1.0:10.0);

julia> y = generalisedlogistic(t; A=0.0, K=10.0, C=1.0, Q=1.0, B=0.75, v=0.1);

julia> (y[1] == minimum(y)) && (y[end] == maximum(y))
true
```
"""
function generalisedlogistic(
    t::Vector{Float64};
    A = 0.0,
    K = 1.0,
    C = 1.0,
    Q = 1.0,
    B = 1.0,
    v = 1.0,
)::Vector{Float64}
    A .+ ((K - A) ./ (C .+ (Q .* exp.(-B .* t))) .^ (1/v))
end

"""
    generalisedlogistic(growth_model::GrowthModel; t::Vector{Float64})::Vector{Float64}

Computes the generalized logistic function for a given vector of time points `t` using the parameters defined in the `GrowthModel`.

# Arguments
- `growth_model::GrowthModel`: A struct containing the parameters for the generalized logistic function.
- `t::Vector{Float64}`: A vector of time points at which to evaluate the function.

# Returns
- `Vector{Float64}`: A vector of values representing the generalized logistic function evaluated at each time point in `t`.

# Returns
- `Vector{Float64}`: A vector of values representing the generalized logistic function evaluated at each time point in `t`.

# Example
```jldoctest; setup = :(using CropGrowth)
julia> t = collect(0.0:1.0:10.0);

julia> growth_model = GrowthModel(A=0.0, K=10.0, C=1.0, Q=1.0, B=0.75, v=0.1);

julia> y = generalisedlogistic(growth_model, t=t);

julia> (y[1] == minimum(y)) && (y[end] == maximum(y))
true
```
"""
function generalisedlogistic(growth_model::GrowthModel; t::Vector{Float64})::Vector{Float64}
    generalisedlogistic(
        t,
        A = growth_model.A,
        K = growth_model.K,
        C = growth_model.C,
        Q = growth_model.Q,
        B = growth_model.B,
        v = growth_model.v,
    )
end

"""
    modelgrowth(; 
        y::Vector{Float64}, 
        t::Vector{Float64}, 
        θ_search_space::Dict{String, Dict{Symbol, Float64}}, 
        maxiters::Int64=10_000, 
        seed::Int64=42, 
        verbose::Bool=false
    )::GrowthModel

Fit a generalized logistic growth model to the given data `y` over time `t`.

# Arguments
- `y::Vector{Float64}`: The observed data points representing the growth values.
- `t::Vector{Float64}`: The corresponding time points for the observed data.
- `θ_search_space::Dict{String, Dict{Symbol, Float64}}`: A dictionary defining the search space for each parameter of the generalized logistic model. Each key corresponds to a parameter name (`"A"`, `"K"`, `"C"`, `"Q"`, `"B"`, `"v"`), and the value is another dictionary with keys `:init`, `:lower`, and `:upper` specifying the initial value, lower bound, and upper bound for that parameter.
- `maxiters::Int64=10_000`: The maximum number of iterations for the optimization algorithm. Defaults to 10,000.
- `seed::Int64=42`: The random seed for reproducibility. Defaults to 42.
- `verbose::Bool=false`: If `true`, prints the fitted parameters, fit statistics, and displays a plot of the fitted model. Defaults to `false`.

# Returns
- `GrowthModel`: A structure containing the fitted parameters of the generalized logistic growth model and fit statistics.

# Details
- The function uses an optimization algorithm to minimize the mean squared error between the observed data `y` and the generalized logistic model. 
- The model parameters are constrained within bounds specified in the `θ_search_space` dictionary.
- The optimization is performed using the `BBO_adaptive_de_rand_1_bin_radiuslimited()` algorithm.
- The function also computes fit statistics, which are included in the returned `GrowthModel` structure.
- If `verbose` is set to `true`, the function prints the fitted parameters, fit statistics, and displays a scatter plot of the observed data along with the fitted curve.

# Notes
- `y_max` may be unexpectedly high if `C` is fitted freely which may indicate a simple linear growth rather than a logistic growth pattern.
- If `y_max` is desired to be equal to `K`, set `C = 1.0` when defining the `θ_search_space`, i.e. In `θ_search_space` define `"C" => Dict(:init=>1.0, :lower=>1.0, :upper=>1.0)`.

# Example
```jldoctest; setup = :(using CropGrowth)
julia> t = collect(1.0:10.0);

julia> growth_model_1 = modelgrowth(t=t, y=sort(vcat([0.0, 0.0, 1.0, 1.0], rand(6))));

julia> growth_model_2 = modelgrowth(t=t, y=rand(10));

julia> growth_model_1.fit_statistics["R²"] > growth_model_2.fit_statistics["R²"]
true
```
"""
function modelgrowth(;
    t::Vector{Float64},
    y::Vector{Float64},
    θ_search_space::Dict{String,Dict{Symbol,Float64}} = Dict(
        "A" => Dict(:init=>minimum(y), :lower=>0.0, :upper=>maximum(y)),
        "K" => Dict(:init=>maximum(y), :lower=>0.0, :upper=>2*maximum(y)),
        "C" => Dict(:init=>1.0, :lower=>0.0, :upper=>5.0),
        "Q" => Dict(:init=>1.0, :lower=>0.0, :upper=>5.0),
        "B" => Dict(:init=>1.0, :lower=>0.0, :upper=>10.0), # assumes no negative growths
        "v" => Dict(:init=>1.0, :lower=>1e-5, :upper=>10.0),
    ),
    maxiters::Int64 = 10_000,
    seed::Int64 = 42,
    verbose::Bool = false,
)::GrowthModel
    # y = sort(vcat([0.0, 0.0, 1.0, 1.0], rand(6))); t = collect(1.0:10.0); maxiters=10_000; seed=42; verbose=false
    # θ_search_space = Dict(
    #     "A" => Dict(:init=>minimum(y), :lower=>0.0, :upper=>maximum(y)),
    #     "K" => Dict(:init=>maximum(y), :lower=>0.0, :upper=>2*maximum(y)),
    #     "C" => Dict(:init=>1.0, :lower=>0.0, :upper=>5.0),
    #     "Q" => Dict(:init=>1.0, :lower=>0.0, :upper=>5.0),
    #     "B" => Dict(:init=>1.0, :lower=>0.0, :upper=>10.0), # assumes no negative growths
    #     "v" => Dict(:init=>1.0, :lower=>1e-5, :upper=>10.0),
    # )
    if !all(haskey(θ_search_space, param) for param in ["A", "K", "C", "Q", "B", "v"])
        throw(
            ArgumentError(
                "θ_search_space must contain keys for all parameters: \"A\", \"K\", \"C\", \"Q\", \"B\", \"v\".",
            ),
        )
    end
    for (k, v) in θ_search_space
        if !all(haskey(v, s) for s in [:init, :lower, :upper])
            throw(
                ArgumentError(
                    "Each parameter in θ_search_space must have :init, :lower, and :upper keys.",
                ),
            )
        end
    end
    opt_fun(θ, t) = mean(
        (
            y .- generalisedlogistic(
                t,
                A = θ[1],
                K = θ[2],
                C = θ[3],
                Q = θ[4],
                B = θ[5],
                v = θ[6],
            )
        ) .^ 2,
    )
    prob = OptimizationProblem(
        opt_fun,
        [
            θ_search_space["A"][:init],
            θ_search_space["K"][:init],
            θ_search_space["C"][:init],
            θ_search_space["Q"][:init],
            θ_search_space["B"][:init],
            θ_search_space["v"][:init],
        ],
        t,
        lb = [
            θ_search_space["A"][:lower],
            θ_search_space["K"][:lower],
            θ_search_space["C"][:lower],
            θ_search_space["Q"][:lower],
            θ_search_space["B"][:lower],
            θ_search_space["v"][:lower],
        ],
        ub = [
            θ_search_space["A"][:upper],
            θ_search_space["K"][:upper],
            θ_search_space["C"][:upper],
            θ_search_space["Q"][:upper],
            θ_search_space["B"][:upper],
            θ_search_space["v"][:upper],
        ],
    )
    θ = solve(
        prob,
        BBO_adaptive_de_rand_1_bin_radiuslimited(),
        maxiters = maxiters,
        seed = seed,
    )
    growth_model = GrowthModel(A = θ[1], K = θ[2], C = θ[3], Q = θ[4], B = θ[5], v = θ[6])
    growth_model.fit_statistics = fitstatistics(
        y = y,
        ŷ = generalisedlogistic(
            t,
            A = growth_model.A,
            K = growth_model.K,
            C = growth_model.C,
            Q = growth_model.Q,
            B = growth_model.B,
            v = growth_model.v,
        ),
    )
    if verbose
        println("Fitted parameters:")
        println("A = $(growth_model.A)")
        println("K = $(growth_model.K)")
        println("C = $(growth_model.C)")
        println("Q = $(growth_model.Q)")
        println("B = $(growth_model.B)")
        println("v = $(growth_model.v)")
        println("y_max = $(growth_model.y_max)")
        println("Fit statistics:")
        for (stat_name, stat_value) in growth_model.fit_statistics
            println("$stat_name = $stat_value")
        end
        t̂ = collect(minimum(t):0.1:maximum(t))
        ŷ = generalisedlogistic(growth_model, t = t̂)
        p = UnicodePlots.scatterplot(t, y, marker = :diamond)
        UnicodePlots.lineplot!(p, t̂, ŷ)
        display(p)
    end
    return growth_model
end

"""
    timetomaxperc(growth_model::GrowthModel; p::Vector{Float64} = [0.5])::Vector{Float64}

Calculate the time required for a growth model to reach a proportion `p` of its carrying capacity `K`.

# Arguments
- `growth_model::GrowthModel`: An instance of the `GrowthModel` type containing the parameters of the growth model.
- `p::Vector{Float64}`: A vector of proportions (default is `[0.5]`) representing the fraction of the carrying capacity `K` to compute the time for.

# Returns
- `Vector{Float64}`: A vector of times corresponding to each proportion in `p`.

# Details
- The function computes the maximum value of the growth model (``y_{max} = A + (K-A)/C^(1/v)``) based on its parameters.
- For each proportion in ``p``, the corresponding value of ``y`` is calculated as ``p * y_{max}``.
- The time required to reach each ``y`` is computed using the generalised logistic growth model formula, which involves logarithmic and complex arithmetic operations.
- The result is the real part of the computed times.
- The formula is:
```math
t = -\\frac{1}{B} \\cdot \\ln \\left( \\frac{\\left( \\left( \\frac{K - A}{y - A} \\right)^v - C \\right)}{Q} \\right)
```

# Notes
- The function assumes that the `GrowthModel` instance contains the parameters `A`, `K`, `C`, `Q`, `B`, and `v`.
- The computation may involve complex numbers, but only the real part of the result is returned.

# Example
```jldoctest; setup = :(using CropGrowth)
julia> growth_model = modelgrowth(t=collect(1.0:10.0), y=sort(vcat([0.0, 0.0, 1.0, 1.0], rand(6))));

julia> p = [0.5, 0.95, 0.99];

julia> time_to_Kp = timetomaxperc(growth_model; p=p);

julia> ŷ = generalisedlogistic(growth_model, t=time_to_Kp);

julia> all( abs.(ŷ .- (p .* growth_model.y_max)) .< 1e-5 )
true
```
"""
function timetomaxperc(
    growth_model::GrowthModel;
    p::Vector{Float64} = [0.5],
)::Vector{Float64}
    # growth_model = modelgrowth(y=sort(vcat([0.0, 0.0, 1.0, 1.0], rand(6))), t=collect(1.0:10.0), maxiters=10_000, seed=42, verbose=true); p = [0.5, 0.9]
    y_max =
        growth_model.A + (growth_model.K-growth_model.A)/growth_model.C^(1/growth_model.v)
    y = p .* y_max
    α = -(1.00/growth_model.B)
    β = (growth_model.K - growth_model.A) ./ (y .- growth_model.A)
    γ = Complex.(β) .^ growth_model.v
    δ = (γ .- growth_model.C) ./ growth_model.Q
    time_to_Kp = [x.re for x in α .* log.(δ)]
    time_to_Kp
end
