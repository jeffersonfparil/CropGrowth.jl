using WGLMakie, ColorSchemes
using Bonito, Markdown
using DataFrames, CSV, StatsBase, Distributions, Random
using WordLists

function simulate()
    Random.seed!(6767)
    years = string.("Year_", 2026:2031)
    seasons = ["Early_Spring", "Late_Spring", "Summer", "Autumn", "Winter"]
    harvests = ["Harvest_1", "Harvest_2", "Harvest_3"]
    sites = lowercase.(rand(words("Dutch"), 5))
    treatments = lowercase.(rand(words("German"), 3))
    entries = lowercase.(rand(words("Spanish"), 25))
    effects_year = rand(Chisq(5), length(years))
    effects_season = [4, 5, 3, 2, 1]
    effects_harvest = rand(Chisq(5), length(harvests))
    effects_site = rand(Chisq(5), length(sites))
    effects_treatment = rand(Chisq(5), length(treatments))
    effects_entry = rand(Chisq(5), length(entries))
    g_x_e = String[]
    yield = Float64[]
    for (i, year) in enumerate(years)
        for (j, season) in enumerate(seasons)
            for (k, harvest) in enumerate(harvests)
                (harvest != harvests[1]) && (rand() < 0.5) ? continue : nothing
                for (l, site) in enumerate(sites)
                    for (m, treatment) in enumerate(treatments)
                        for (n, entry) in enumerate(entries)
                            # i=1; j=1; k=1; l=1; m=1; n=1; year = years[i]; season = seasons[j]; harvest = harvests[k]; site = sites[l]; treatment = treatments[m]; entry = entries[n];
                            push!(g_x_e, join([year, season, harvest, site, treatment, entry], "\t"))
                            push!(
                                yield,
                                effects_year[i] + 
                                effects_season[j] + 
                                effects_harvest[k] + 
                                effects_site[l] + 
                                effects_treatment[m] + 
                                effects_entry[n] + 
                                rand(Normal(0, 1))
                            )
                        end
                    end
                end
            end
        end
    end
    df = DataFrame(
        years = [String(x[1]) for x in split.(g_x_e, "\t")],
        seasons = [String(x[2]) for x in split.(g_x_e, "\t")],
        harvests = [String(x[3]) for x in split.(g_x_e, "\t")],
        sites = [String(x[4]) for x in split.(g_x_e, "\t")],
        treatments = [String(x[5]) for x in split.(g_x_e, "\t")],
        entries = [String(x[6]) for x in split.(g_x_e, "\t")],
        yield = 5 .* yield,
    )
    df
end

df = simulate()

WGLMakie.activate!(; use_html_widgets = true)

fig = Figure(size=(1_500, 900))
layout = GridLayout(fig[1, 1], tellwidth = false)

year = unique(df.years)[1]
site = unique(df.sites)[1]
treatment = unique(df.treatments)[1]
df_tmp = df |>
    d -> filter(x -> (x.years == year) && (x.sites == site) && (x.treatments == treatment), d) |> 
    d -> groupby(d, [:years, :seasons, :entries]) |>
    d -> combine(d, "yield" => sum => "yield")
df_tmp.x_entries = [findfirst(unique(df_tmp.entries) .== x) for x in df_tmp.entries]
df_tmp.x_seasons = [findfirst(unique(df_tmp.seasons) .== x) for x in df_tmp.seasons]
n = length(unique(df_tmp.entries))
m = length(unique(df_tmp.seasons))
colours_entries = resample(ColorSchemes.tol_bright, n)
colours_seasons = resample(ColorSchemes.tol_sunset, m)

axis_stacked_barplot = Axis(
    layout[1,1],
    ylabel = "Yield",
    xticks = (1:n, unique(df_tmp.entries)),
    xticklabelrotation = π/(1 + (75/n)),
)
barplot!(
    axis_stacked_barplot,
    df_tmp.x_entries,
    df_tmp.yield,
    stack = df_tmp.x_seasons,
    color = colours_seasons[df_tmp.x_seasons],
)

axis_dodged_barplot = Axis(
    layout[1,2],
    ylabel = "Yield",
    xticks = (1:n, unique(df_tmp.entries)),
    xticklabelrotation = π/(1 + (75/n)),
)
barplot!(
    axis_dodged_barplot,
    df_tmp.x_entries,
    df_tmp.yield,
    dodge = df_tmp.x_seasons,
    color = colours_seasons[df_tmp.x_seasons],
)

axes = Dict()
for (i, season) in enumerate(unique(df_tmp.seasons))
    # i = 1; season = unique(df_tmp.seasons)[i]
    row = Int(ceil((i+2) / 2))
    col = ((i+1) % 2) + 1
    # println("row: $row; col: $col")
    axes[i] = Axis(
        layout[row, col],
        title = season,
        ylabel = "Yield",
        xticks = (1:n, unique(df_tmp.entries)),
        xticklabelrotation = π/(1 + (75/n)),
    )
    idx = findall(df_tmp.seasons .== season)
    barplot!(
        axes[i],
        df_tmp.x_entries[idx],
        df_tmp.yield[idx],
        color = colours_seasons[i],
    )
end

Legend(
    layout[:,0],
    [[MarkerElement(color=x, marker=:circle, markersize=15, strokecolor=:black)] for x in colours_seasons],
    unique(df_tmp.seasons),
)

