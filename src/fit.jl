"""
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

Fits generalized logistic growth models to the data provided in the input `DataFrame` and returns a `DataFrame` containing the fitted parameters, fit statistics, and time to reach specified percentages of the final value.

# Arguments
- `df::DataFrame`: Input data containing the required columns specified in `REQUIRED_COLUMNS` and at least one trait column.
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

# Returns
- `DataFrame`: A `DataFrame` containing the fitted parameters (`A`, `K`, `C`, `Q`, `B`, `v`), fit statistics, value of the growth models at ``t=0`` (`y_t0`), maximum value of the growth model (`y_max`), and time to reach specified percentages of the final value for each combination of entry, site, replication, and growing period.

# Notes
- The input `DataFrame` must contain the required columns specified in the global variable `REQUIRED_COLUMNS`, as well as at least one additional trait column.
- If the `DataFrame` contains more than one trait column, only the first trait column will be used.
- Combinations with fewer than `min_t` time points will be skipped, and a warning will be issued.
- The function uses a progress bar to indicate the fitting process if `verbose=true`.

# Details
Fits a generalized logistic growth model to the data using the following equation:

```math
y(t) = A + \\frac{K-A}{C + (Qe^{-Bt})^{1/v}}
```

where:

- ``y(t)``: biomass at time ``t``
- ``A``: lower asymptote (initial or minimum biomass)
- ``K``: positively affects the upper asymptote (can be the final or maximum biomass if ``C = 1.00``, since ``y_{max} = A + (K-A)/C^(1/v)``, then ``y_{max} = A + K - A``, therefore: ``y_{max} = K``)
- ``C``: negatively affects the final or maximum biomass
- ``Q``: negatively affects initial or minimum biomass
- ``e``: Euler's number (~2.71828)
- ``B``: growth rate
- ``v``: asymmetry parameter (``v ≥ 0``; small values: fast growth early; large values: fast growth later)

Additional information are provided in the output `DataFrame`:
- `y_t0` is the value of the growth model at time ``t = 0``
- `y_max` is the maximum value of the growth model (``y_{max} = A + (K-A)/C^(1/v)``)
- One of the following fit statistic such as R² (default), RMSE, MSE, MAE, and Pearson's correlation coefficient (ρ)

# Example
```jldoctest; setup = :(using CropGrowth, DataFrames, StatsBase)
julia> df = simulate(n_entries=5, seed=42);

julia> df_out, skipped_combinations = fitgrowthmodels(df);

julia> length(unique(df.entries)) == length(unique(df_out.entries))
true

julia> length(unique(df.sites)) == length(unique(df_out.sites))
true

julia> length(unique(df.replications)) == length(unique(df_out.replications))
true

julia> length(unique(df.growing_periods)) == length(unique(df_out.growing_periods))
true

julia> all(x -> x ∈ names(df_out), ["A", "K", "C", "Q", "B", "v", "y_max", "R²", "time_to_50p", "time_to_90p"])
true

julia> (var(df_out.A) > 0.0) && (var(df_out.K) > 0.0) && (var(df_out.C) == 0.0) && (var(df_out.Q) == 0.0) && (var(df_out.B) > 0.0) && (var(df_out.v) > 0.0)
true

julia> (var(df_out.y_max) > 0.0) && (var(df_out."R²") > 0.0) && (var(df_out.time_to_50p) > 0.0) && (var(df_out.time_to_90p) > 0.0)
true
```
"""
function fitgrowthmodels(
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
)::Tuple{DataFrame,Vector{String}}
    # df::DataFrame = simulate(); min_t::Int64=3; perc_of_final::Vector{Float64} = [0.5, 0.9]; fit_statistic::String = "R²"; maxiters::Int64 = 10_000; seed::Int64 = 42; show_plots::Bool=false; verbose::Bool=false
    # A = Dict(:init=>minimum(select(df, Not(REQUIRED_COLUMNS))[:, 1]), :lower=>0.0, :upper=>maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1])); K = Dict(:init=>maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]), :lower=>0.0, :upper=>2*maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1])); C = Dict(:init=>1.0, :lower=>1.0, :upper=>1.0); Q = Dict(:init=>1.0, :lower=>1.0, :upper=>1.0); B = Dict(:init=>1.0, :lower=>0.0, :upper=>10.0); v = Dict(:init=>1.0, :lower=>1e-5, :upper=>10.0)
    if !all(x -> x ∈ names(df), REQUIRED_COLUMNS)
        throw(
            ArgumentError(
                "DataFrame is missing required columns. Required columns are: $(REQUIRED_COLUMNS)",
            ),
        )
    end
    if ncol(df) < length(REQUIRED_COLUMNS) + 1
        throw(
            ArgumentError(
                "DataFrame must contain at least one trait column in addition to the required columns: $(REQUIRED_COLUMNS)",
            ),
        )
    end
    if ncol(df) > length(REQUIRED_COLUMNS) + 1
        @warn "DataFrame contains more than one trait column. Only the first trait column will be used, i.e. \"$(setdiff(names(df), REQUIRED_COLUMNS)[1])\"."
    end
    for p in [A, K, C, Q, B, v]
        if !all(haskey(p, s) for s in [:init, :lower, :upper])
            throw(
                ArgumentError(
                    "Parameter search space for $p must have :init, :lower, and :upper keys.",
                ),
            )
        end
    end
    if fit_statistic ∉ ["R²", "RMSE", "MSE", "MAE", "ρ"]
        throw(
            ArgumentError(
                "fit_statistic must be one of: [\"R²\", \"RMSE\", \"MSE\", \"MAE\", \"ρ\"]",
            ),
        )
    end
    θ_search_space::Dict{String,Dict{Symbol,Float64}} =
        Dict("A" => A, "K" => K, "C" => C, "Q" => Q, "B" => B, "v" => v)
    entries = sort(unique(df[!, :entries]))
    sites = sort(unique(df[!, :sites]))
    replications = sort(unique(df[!, :replications]))
    growing_periods = sort(unique(df[!, :growing_periods]))
    trait_name = setdiff(names(df), REQUIRED_COLUMNS)[1]
    skipped_combinations::Vector{String} = []
    fitted_parameters = Dict(
        "entries" => String[],
        "sites" => String[],
        "replications" => String[],
        "growing_periods" => String[],
        "A" => Float64[],
        "K" => Float64[],
        "C" => Float64[],
        "Q" => Float64[],
        "B" => Float64[],
        "v" => Float64[],
        "y_t0" => Float64[],
        "y_max" => Float64[],
        fit_statistic => Float64[],
    )
    for p in perc_of_final
        fitted_parameters["time_to_$(Int(round(p*100)))p"] = Float64[]
    end
    if verbose
        pb = ProgressMeter.Progress(
            length(entries) *
            length(sites) *
            length(replications) *
            length(growing_periods),
            desc = "Fitting growth models",
        )
    end
    for entry in entries
        for site in sites
            for replication in replications
                for growing_period in growing_periods
                    # entry = entries[1]; site = sites[1]; replication = replications[1]; growing_period = growing_periods[1];
                    idx = findall(
                        (df.entries .== entry) .&&
                        (df.sites .== site) .&&
                        (df.replications .== replication) .&&
                        (df.growing_periods .== growing_period),
                    )
                    combination = "entry=\"$entry\"; site=\"$site\"; replication=\"$replication\"; growing_period=\"$growing_period\""
                    if length(idx) < min_t
                        # @warn "Not enough data points (minimum t = $min_t) to fit growth model for $combination. Skipping."
                        push!(skipped_combinations, combination)
                        continue
                    end
                    df_sub = df[idx, :]
                    t::Vector{Float64} = df_sub.time_points
                    y::Vector{Float64} = df_sub[:, trait_name]
                    if show_plots
                        println("Fitting growth model for $combination")
                    end
                    growth_model = modelgrowth(
                        y = y,
                        t = t,
                        θ_search_space = θ_search_space,
                        maxiters = maxiters,
                        seed = seed,
                        verbose = show_plots,
                    )
                    time_to_perc_of_final = timetomaxperc(growth_model, p = perc_of_final)
                    fitted_parameters["entries"] = vcat(fitted_parameters["entries"], entry)
                    fitted_parameters["sites"] = vcat(fitted_parameters["sites"], site)
                    fitted_parameters["replications"] =
                        vcat(fitted_parameters["replications"], replication)
                    fitted_parameters["growing_periods"] =
                        vcat(fitted_parameters["growing_periods"], growing_period)
                    fitted_parameters["A"] = vcat(fitted_parameters["A"], growth_model.A)
                    fitted_parameters["K"] = vcat(fitted_parameters["K"], growth_model.K)
                    fitted_parameters["C"] = vcat(fitted_parameters["C"], growth_model.C)
                    fitted_parameters["Q"] = vcat(fitted_parameters["Q"], growth_model.Q)
                    fitted_parameters["B"] = vcat(fitted_parameters["B"], growth_model.B)
                    fitted_parameters["v"] = vcat(fitted_parameters["v"], growth_model.v)
                    fitted_parameters["y_t0"] =
                        vcat(fitted_parameters["y_t0"], growth_model.y_t0)
                    fitted_parameters["y_max"] =
                        vcat(fitted_parameters["y_max"], growth_model.y_max)
                    fitted_parameters[fit_statistic] = vcat(
                        fitted_parameters[fit_statistic],
                        growth_model.fit_statistics[fit_statistic],
                    )
                    for (i, p) in enumerate(perc_of_final)
                        # i = 1; p = perc_of_final[i]
                        Kp = time_to_perc_of_final[i]
                        fitted_parameters["time_to_$(Int(round(p*100)))p"] =
                            vcat(fitted_parameters["time_to_$(Int(round(p*100)))p"], Kp)
                    end
                    if verbose
                        ProgressMeter.next!(pb)
                    end
                end
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
        println(
            "Skipped $(length(skipped_combinations)) combinations due to insufficient data points (i.e. t < $min_t).",
        )
    end
    # Output dataFrame and sort columns sensibly
    df_out = DataFrame(fitted_parameters)
    select!(
        df_out,
        :entries,
        :sites,
        :replications,
        :growing_periods,
        :A,
        :K,
        :C,
        :Q,
        :B,
        :v,
        :y_t0,
        :y_max,
        Symbol(fit_statistic),
        Not([
            :entries,
            :sites,
            :replications,
            :growing_periods,
            :A,
            :K,
            :C,
            :Q,
            :B,
            :v,
            :y_max,
            Symbol(fit_statistic),
        ]),
    )
    sort!(df_out, [:entries, :sites, :replications, :growing_periods])
    return (df_out, skipped_combinations)
end
