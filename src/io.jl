"""
    simulate(; 
        n_entries::Int64 = 100, 
        n_sites::Int64 = 2, 
        n_replications::Int64 = 3, 
        n_growing_periods::Int64 = 4, 
        n_time_points_per_growing_period::Int64 = 10, 
        seed::Int64 = 42
    )::DataFrame

Simulates crop growth data and returns it as a `DataFrame`. The function generates synthetic data for multiple entries, sites, replications, growing periods, and time points, with random biomass values.

# Keyword Arguments
- `n_entries::Int64`: Number of crop entries to simulate (default: 100).
- `n_sites::Int64`: Number of sites to simulate (default: 2).
- `n_replications::Int64`: Number of replications per site (default: 3).
- `n_growing_periods::Int64`: Number of growing periods per replication (default: 4).
- `n_time_points_per_growing_period::Int64`: Number of time points per growing period (default: 10).
- `seed::Int64`: Random seed for reproducibility (default: 42).

# Returns
A `DataFrame` with the following columns:
- `entries`: Crop entry identifiers.
- `sites`: Site identifiers.
- `replications`: Replication identifiers.
- `growing_periods`: Growing period identifiers.
- `time_points`: Time points within each growing period.
- `biomass`: Simulated biomass values.

# Example
```jldoctest; setup = :(using CropGrowth, DataFrames, StatsBase)
julia> df = simulate(n_entries=5, n_sites=1, n_replications=1, n_growing_periods=1, n_time_points_per_growing_period=5, trait_name="some_trait_name", seed=123);

julia> size(df)
(25, 6)

julia> var(df.some_trait_name) > 0.0
true
```
"""
function simulate(;
    n_entries::Int64 = 100,
    n_sites::Int64 = 2,
    n_replications::Int64 = 3,
    n_growing_periods::Int64 = 4,
    n_time_points_per_growing_period::Int64 = 10,
    trait_name::String = "biomass",
    seed::Int64 = 42,
)::DataFrame
    # n_entries::Int64=100; n_sites::Int64=2; n_replications::Int64=3; n_growing_periods::Int64=4; n_time_points_per_growing_period::Int64=10; trait_name::String = "biomass"; seed::Int64 = 42
    Random.seed!(seed)
    entries = [String("entry_$i") for i = 1:n_entries]
    sites = [String("site_$i") for i = 1:n_sites]
    replications = [String("replication_$i") for i = 1:n_replications]
    growing_periods = [String("growing_period_$i") for i = 1:n_growing_periods]
    time_points = [Float64(i) for i = 0:(n_time_points_per_growing_period-1)]
    out_entries = []
    out_sites = []
    out_replications = []
    out_growing_periods = []
    out_time_points = []
    out_biomass = []
    for entry in entries
        for site in sites
            for replication in replications
                for growing_period in growing_periods
                    for time_point in time_points
                        push!(out_entries, entry)
                        push!(out_sites, site)
                        push!(out_replications, replication)
                        push!(out_growing_periods, growing_period)
                        push!(out_time_points, time_point)
                        push!(
                            out_biomass,
                            maximum([time_point*rand(), time_point*rand()^2]),
                        )
                    end
                end
            end
        end
    end
    df = DataFrame(
        entries = out_entries,
        sites = out_sites,
        replications = out_replications,
        growing_periods = out_growing_periods,
        time_points = out_time_points,
        biomass = out_biomass,
    )
    rename!(df, :biomass => trait_name)
    df
end

"""
    readdelimited(; fname::String, delim::Union{Char,String} = "\t", trait_name::String = "biomass")::DataFrame

Reads a delimited file into a `DataFrame` and validates its structure.

# Arguments
- `fname::String`: The path to the file to be read.
- `delim::Union{Char,String}`: The delimiter used in the file. Defaults to tab (`"\t"`).
- `trait_name::String`: The name of the trait column to validate. Defaults to `"biomass"`.

# Returns
- A `DataFrame` containing the data from the file with the following required columns:
    + `entries`
    + `sites`
    + `replications`
    + `growing_periods`
    + `time_points`
    + The column specified by `trait_name`

# Throws
- `ArgumentError`: If the file does not exist.
- `ArgumentError`: If any of the required columns, including `trait_name`, are missing.

# Notes
- The function expects the file to contain specific required columns, defined in `REQUIRED_COLUMNS`, along with the column specified by `trait_name`.

# Example
```jldoctest; setup = :(using CropGrowth, DataFrames, CSV)
julia> df = simulate();

julia> CSV.write("test.tsv", df, delim="\t");

julia> df_read = readdelimited(fname="test.tsv", delim="\t", trait_name="biomass");

julia> isequal(df, df_read)
true
```
"""
function readdelimited(;
    fname::String,
    delim::Union{Char,String} = "\t",
    trait_name::String = "biomass",
)::DataFrame
    # fname::String = "test.tsv"; df = simulate(); CSV.write(fname, df, delim="\t"); delim::Union{Char, String} = "\t"; trait_name::String="biomass"
    if !isfile(fname)
        throw(ArgumentError("File does not exist: \"$fname\""))
    end
    df = CSV.read(fname, DataFrame, delim = delim)
    expected_column_names = vcat(REQUIRED_COLUMNS, trait_name)
    for expected_column_name in expected_column_names
        if expected_column_name ∉ names(df)
            throw(ArgumentError("Missing column: \"$expected_column_name\""))
        end
    end
    df
end

"""
    writedelimited(df::DataFrame; fname::String, delim::Union{Char,String} = "\t", overwrite::Bool = false)::Nothing

Writes the contents of a `DataFrame` to a delimited text file.

# Arguments
- `df::DataFrame`: The data frame to be written to the file. Must contain the required columns:
    + `entries`
    + `sites`
    + `replications`
    + `growing_periods`
    + `time_points`
    + 1 additional column representing the trait data
- `fname::String`: The name of the output file, including its path if necessary.
- `delim::Union{Char,String}`: The delimiter to use in the output file. Defaults to a tab character (`"\t"`).
- `overwrite::Bool`: Whether to overwrite the file if it already exists. Defaults to `false`.

# Returns
- `Nothing`: This function does not return a value.

# Example
```jldoctest; setup = :(using CropGrowth, DataFrames)
julia> df = simulate();

julia> fname = writedelimited(df, fname="test.tsv", delim="\t", overwrite=true);

julia> isfile(fname)
true

julia> df_read = readdelimited(fname="test.tsv", delim="\t", trait_name="biomass");

julia> isequal(df, df_read)
true
```
"""
function writedelimited(
    df::DataFrame;
    fname::String,
    delim::Union{Char,String} = "\t",
    overwrite::Bool = false,
)::String
    # df::DataFrame = simulate(); fname::String="test.tsv"; delim::Union{Char, String}="\t"; overwrite::Bool = false
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
    if isfile(fname) && !overwrite
        throw(
            ArgumentError(
                "File already exists: \"$fname\". Set `overwrite=true` to overwrite.",
            ),
        )
    end
    CSV.write(fname, df, delim = delim)
    fname
end
