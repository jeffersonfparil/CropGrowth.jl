function simulate(;
    n_entries::Int64=100,
    n_sites::Int64=2,
    n_replications::Int64=3,
    n_growing_periods::Int64=4,
    n_time_points_per_growing_period::Int64=10,
)::DataFrame
    # n_entries::Int64=100; n_sites::Int64=2; n_replications::Int64=3; n_growing_periods::Int64=4; n_time_points_per_growing_period::Int64=10
    entries = [String("entry_$i") for i in 1:n_entries]
    sites = [String("site_$i") for i in 1:n_sites]
    replications = [String("replication_$i") for i in 1:n_replications]
    growing_periods = [String("growing_period_$i") for i in 1:n_growing_periods]
    time_points = [Float64(i) for i in 0:(n_time_points_per_growing_period-1)]
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
                        push!(out_biomass, time_point*(rand()^2))
                    end
                end
            end
        end
    end
    DataFrame(
        entries = out_entries,
        sites = out_sites,
        replications = out_replications,
        growing_periods = out_growing_periods,
        time_points = out_time_points,
        biomass = out_biomass,
    )
end

function readdelimited(; fname::String, delim::Union{Char, String}="\t", trait_name::String="biomass")::DataFrame
    # fname::String = begin
    #     df = simulate()
    #     fname = "test.tsv"
    #     CSV.write(fname, df, delim="\t")
    # end
    # delim::Union{Char, String} = "\t"
    # trait_name::String="biomassa"
    df = CSV.read(fname, DataFrame, delim=delim)
    expected_column_names = ["entries", "sites", "replications", "growing_periods", "time_points", trait_name]
    for expected_column_name in expected_column_names
        if expected_column_name ∉ names(df)
            throw(ArgumentError("Missing column: \"$expected_column_name\""))
        end
    end
    df    
end