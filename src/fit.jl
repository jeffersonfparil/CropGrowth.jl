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
        frac_of_final::Vector{Float64} = [0.5, 0.9],
        fit_statistic::String = "R²",
        maxiters::Int64 = 10_000,
        seed::Int64 = 42,
        show_plots::Bool = false,
        verbose::Bool = false,
    )::Tuple{DataFrame, Vector{String}}

Fits generalised logistic growth models to the data provided in the input `DataFrame` and returns a `DataFrame` containing the fitted parameters, fit statistics, and time to reach specified fractions of the final value.

# Arguments
- `df::DataFrame`: Input data containing the required columns specified in `REQUIRED_COLUMNS = ["entries", "sites", "replications", "growing_periods", "time_points"]` and at least one trait column.
- `A::Dict`: Search space for the parameter `A` (lower asymptote). Contains `:init`, `:lower`, and `:upper` keys. Defaults to the minimum and maximum of the trait column with `init=minimum`.
- `K::Dict`: Search space for the parameter `K` (upper asymptote). Contains `:init`, `:lower`, and `:upper` keys. Defaults to the minimum and 2×maximum of the trait column with `init=maximum`.
- `C::Dict`: Search space for the parameter `C`. Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=1.0`, `upper=1.0`.
- `Q::Dict`: Search space for the parameter `Q`. Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=1.0`, `upper=1.0`.
- `B::Dict`: Search space for the parameter `B` (growth rate). Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=0.0`, `upper=10.0`.
- `v::Dict`: Search space for the parameter `v` (asymmetry parameter). Contains `:init`, `:lower`, and `:upper` keys. Defaults to `init=1.0`, `lower=1e-5`, `upper=10.0`.
- `min_t::Int64`: Minimum number of time points required to fit the growth model for a specific combination of entry, site, replication, and growing period. Defaults to `3`.
- `frac_of_final::Vector{Float64}`: Percentages of the final value for which the time to reach these fractions will be calculated. Defaults to `[0.5, 0.9]`.
- `fit_statistic::String`: The fit statistic to be used for evaluating the model. Must be one of `["R²", "RMSE", "MSE", "MAE", "ρ"]`. Defaults to `"R²"`.
- `maxiters::Int64`: Maximum number of iterations allowed for the optimisation process. Defaults to `10_000`.
- `seed::Int64`: Random seed for reproducibility. Defaults to `42`.
- `show_plots::Bool`: Whether to show fitted growth curve plots. Defaults to `false`.
- `verbose::Bool`: Whether to display progress and additional information during the fitting process. Defaults to `false`.

# Returns
- `Tuple{DataFrame, Vector{String}}`: 
    - The first element is a `DataFrame` containing the fitted parameters (`A`, `K`, `C`, `Q`, `B`, `v`), fit statistics, value of the growth models at ``t=0`` (`y_t0`), maximum value of the growth model (`y_max`), time to reach specified fractions of the final value (`time_to_*p`), and number of time-points (`number_of_time_points`) for each combination of entry, site, replication, and growing period.
    - The second element is a `Vector{String}` containing the combinations that were skipped due to insufficient data points.

# Notes
- The input `DataFrame` must contain the required columns specified in the global variable `REQUIRED_COLUMNS = ["entries", "sites", "replications", "growing_periods", "time_points"]`, as well as at least one additional trait column.
- If the `DataFrame` contains more than one trait column, only the first trait column will be used.
- Combinations with fewer than `min_t` time points will be skipped.
- The function uses a progress bar to indicate the fitting process if `verbose=true`.
- The optimisation is performed using the `BBO_adaptive_de_rand_1_bin_radiuslimited()` algorithm ([details of the optimisation algorithm](https://docs.sciml.ai/Optimization/stable/optimisation_packages/blackboxoptim/)).
- The optimisation algorithm minimises the mean squared error between the observed data `y` and the generalised logistic model.
- Parallel curve fitting:
    - The per-combination fits are executed in parallel using `Threads.@threads` over the prepared list of filter combinations. To take advantage of this, start Julia with multiple threads (for example `julia -t 4` or set `JULIA_NUM_THREADS`).
    - Access to the shared `fitted_parameters` dictionary is protected with a `ReentrantLock` and `@lock` to ensure thread-safety when appending results. This avoids race conditions but introduces some serialization for the final writes.

# Details
Fits a generalised logistic growth model to the data using the following equation:

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
        :init => minimum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
        :lower => 0.0,
        :upper => maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
    ),
    K = Dict(
        :init => maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
        :lower => 0.0,
        :upper => 2 * maximum(select(df, Not(REQUIRED_COLUMNS))[:, 1]),
    ),
    C = Dict(:init => 1.0, :lower => 1.0, :upper => 1.0),
    Q = Dict(:init => 1.0, :lower => 1.0, :upper => 1.0),
    B = Dict(:init => 1.0, :lower => 0.0, :upper => 10.0),
    v = Dict(:init => 1.0, :lower => 1e-5, :upper => 10.0),
    min_t::Int64 = 3,
    frac_of_final::Vector{Float64} = [0.5, 0.9],
    fit_statistic::String = "R²",
    maxiters::Int64 = 10_000,
    seed::Int64 = 42,
    show_plots::Bool = false,
    verbose::Bool = false,
)::Tuple{DataFrame,Vector{String}}
    # df::DataFrame = simulate(); min_t::Int64=3; frac_of_final::Vector{Float64} = [0.5, 0.9]; fit_statistic::String = "R²"; maxiters::Int64 = 10_000; seed::Int64 = 42; show_plots::Bool=false; verbose::Bool=false
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
    # Identify growth curves we can fit, i.e. >= min_t
    curve_ids::Vector{Dict{Symbol,String}} = []
    skipped_combinations::Vector{String} = []
    dict_entry = Dict()
    for entry in entries
        dict_entry[String(entry)] = df.entries .== entry
    end
    dict_site = Dict()
    for site in sites
        dict_site[String(site)] = df.sites .== site
    end
    dict_replication = Dict()
    for replication in replications
        dict_replication[String(replication)] = df.replications .== replication
    end
    dict_growing_period = Dict()
    for growing_period in growing_periods
        dict_growing_period[String(growing_period)] = df.growing_periods .== growing_period
    end
    if verbose
        pb = ProgressMeter.Progress(
            length(dict_entry) *
            length(dict_site) *
            length(dict_replication) *
            length(dict_growing_period),
            desc = "Preparing growth data combinations",
        )
    end
    for (entry, idx_entry) in sort(dict_entry)
        for (site, idx_site) in sort(dict_site)
            for (replication, idx_replication) in sort(dict_replication)
                for (growing_period, idx_growing_period) in sort(dict_growing_period)
                    # entry = string.(keys(dict_entry))[1]; idx_entry = dict_entry[entry]; site = string.(keys(dict_site))[1]; idx_site = dict_site[site]; replication = string.(keys(dict_replication))[1]; idx_replication = dict_replication[replication]; growing_period = string.(keys(dict_growing_period))[1]; idx_growing_period = dict_growing_period[growing_period]
                    idx = idx_entry .&& idx_site .&& idx_replication .&& idx_growing_period
                    if verbose
                        ProgressMeter.next!(pb)
                    end
                    combination = "entry=\"$entry\"; site=\"$site\"; replication=\"$replication\"; growing_period=\"$growing_period\""
                    if sum(idx) < min_t
                        # @warn "Not enough data points (minimum t = $min_t) to fit growth model for $combination. Skipping."
                        push!(skipped_combinations, combination)
                        continue
                    end
                    push!(
                        curve_ids,
                        Dict(
                            :entry => entry,
                            :site => site,
                            :replication => replication,
                            :growing_period => growing_period,
                        ),
                    )
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
    # Prepare output data dictionary
    fitted_parameters = Dict(
        "entries" => String[],
        "sites" => String[],
        "replications" => String[],
        "growing_periods" => String[],
        "number_of_time_points" => Int64[],
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
    for p in frac_of_final
        fitted_parameters["time_to_$(Int(round(p*100)))p"] = Float64[]
    end
    # Fit growth curves in parallel (remember open julia with something like the following to enable multiple threads: `julia +1.12 --threads=23,1 --project=.`)
    n = length(curve_ids)
    if verbose
        pb = ProgressMeter.Progress(n; desc = "Fitting growth models")
    end
    thread_lock::ReentrantLock = ReentrantLock()
    Threads.@threads for i = 1:n
        # i = 1
        entry = curve_ids[i][:entry]
        site = curve_ids[i][:site]
        replication = curve_ids[i][:replication]
        growing_period = curve_ids[i][:growing_period]
        idx = findall(
            (df.entries .== entry) .&&
            (df.sites .== site) .&&
            (df.replications .== replication) .&&
            (df.growing_periods .== growing_period),
        )
        combination = "entry=\"$entry\"; site=\"$site\"; replication=\"$replication\"; growing_period=\"$growing_period\""
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
        time_to_frac_of_final = timetomaxperc(growth_model, p = frac_of_final)
        # Collect results locally to minimise lock contention
        local_entry = entry
        local_site = site
        local_replication = replication
        local_growing_period = growing_period
        local_n_t = length(idx)
        local_A = growth_model.A
        local_K = growth_model.K
        local_C = growth_model.C
        local_Q = growth_model.Q
        local_B = growth_model.B
        local_v = growth_model.v
        local_y_t0 = growth_model.y_t0
        local_y_max = growth_model.y_max
        local_fit_stat = growth_model.fit_statistics[fit_statistic]
        local_times = [time_to_frac_of_final[j] for j = 1:length(frac_of_final)]
        # Append all results with single lock
        @lock thread_lock begin
            push!(fitted_parameters["entries"], local_entry)
            push!(fitted_parameters["sites"], local_site)
            push!(fitted_parameters["replications"], local_replication)
            push!(fitted_parameters["growing_periods"], local_growing_period)
            push!(fitted_parameters["number_of_time_points"], local_n_t)
            push!(fitted_parameters["A"], local_A)
            push!(fitted_parameters["K"], local_K)
            push!(fitted_parameters["C"], local_C)
            push!(fitted_parameters["Q"], local_Q)
            push!(fitted_parameters["B"], local_B)
            push!(fitted_parameters["v"], local_v)
            push!(fitted_parameters["y_t0"], local_y_t0)
            push!(fitted_parameters["y_max"], local_y_max)
            push!(fitted_parameters[fit_statistic], local_fit_stat)
            for (j, p) in enumerate(frac_of_final)
                key = "time_to_$(Int(round(p*100)))p"
                push!(fitted_parameters[key], local_times[j])
            end
            if verbose
                ProgressMeter.next!(pb)
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
    end
    # Output dataFrame and sort columns sensibly
    df_out = DataFrame(fitted_parameters)
    select!(
        df_out,
        :entries,
        :sites,
        :replications,
        :growing_periods,
        :number_of_time_points,
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
            :number_of_time_points,
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
