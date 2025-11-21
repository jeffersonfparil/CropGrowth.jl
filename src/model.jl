"""

    struct GrowthModel

Holds parameters for the generalized logistic growth model:

```math
y(t) = A + \\frac{K-A}{C + (Qe^{-Bt})^{1/v}}
```

where:

- ``y(t)``: biomass at time ``t``
- ``A``: lower asymptote (initial or minimum biomass)
- ``K``: upper asymptote (can be the final or maximum biomass if ``C = 1.00``, since ``y_max = A + (K-A)/C^(1/v)``, then ``y_max = A + K - A``, therefore: ``y_max = K``)
- ``C``: negatively affects the final or maximum biomass
- ``Q``: negatively affects initial or minimum biomass
- ``e``: Euler's number (~2.71828)
- ``B``: growth rate
- ``v``: asymmetry parameter (``v ≥ 0``; small values: fast growth early; large values: fast growth later)
"""
struct GrowthModel
    A::Float64
    K::Float64
    C::Float64
    Q::Float64
    B::Float64
    v::Float64
end

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

function modelgrowth(;
    y::Vector{Float64},
    t::Vector{Float64},
    maxiters::Int64 = 10_000,
    seed::Int64 = 42,
    verbose::Bool = false,
)::GrowthModel
    # y = sort(vcat([0.0, 0.0, 1.0, 1.0], rand(6))); t = collect(1.0:10.0); maxiters=10_000; seed=42; verbose=false
    θ_search_space = Dict(
        "A" => Dict(:init=>minimum(y), :lower=>0.0, :upper=>maximum(y)),
        "K" => Dict(:init=>maximum(y), :lower=>0.0, :upper=>2*maximum(y)),
        "C" => Dict(:init=>1.0, :lower=>0.0, :upper=>5.0),
        "Q" => Dict(:init=>1.0, :lower=>0.0, :upper=>5.0),
        "B" => Dict(:init=>1.0, :lower=>0.0, :upper=>10.0), # assumes no negative growths
        "v" => Dict(:init=>1.0, :lower=>1e-5, :upper=>10.0),
    )
    opt_fun(θ, t) = mean(
        (
            y - generalisedlogistic(
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
    if verbose
        println("Fitted parameters:")
        println("A = $(θ[1])")
        println("K = $(θ[2])")
        println("C = $(θ[3])")
        println("Q = $(θ[4])")
        println("B = $(θ[5])")
        println("v = $(θ[6])")
        t̂ = collect(minimum(t):0.1:maximum(t))
        ŷ = generalisedlogistic(
            t̂,
            A = θ[1],
            K = θ[2],
            C = θ[3],
            Q = θ[4],
            B = θ[5],
            v = θ[6],
        )
        p = UnicodePlots.scatterplot(t, y, marker = :diamond)
        UnicodePlots.lineplot!(p, t̂, ŷ)
        display(p)
    end
    return GrowthModel(θ[1], θ[2], θ[3], θ[4], θ[5], θ[6])
end

function timetoKp(growth_model::GrowthModel; p::Vector{Float64} = [0.5])::Vector{Float64}
    # growth_model = modelgrowth(y=sort(vcat([0.0, 0.0, 1.0, 1.0], rand(6))), t=collect(1.0:10.0), maxiters=10_000, seed=42, verbose=true); p = [0.5, 0.9]
    t =
        -(1.00/growth_model.B) .* log.(
            (
                (
                    Complex.(
                        (growth_model.K - growth_model.A) ./
                        ((p .* growth_model.K) .- growth_model.A),
                    ) .^ growth_model.v
                ) .- growth_model.C
            ) ./ growth_model.Q,
        )
    time_to_Kp = [x.re for x in t]
    @assert sum(
        abs.(generalisedlogistic(growth_model, t = time_to_Kp) .- (p .* growth_model.K)),
    ) < 1e-7
    time_to_Kp
end
