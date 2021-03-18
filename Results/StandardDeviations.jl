#=

    StandardDeviations.jl

    Compute standard deviation for each summary statistic as a function of mean.

=#

using PoresIdentifiability
using StatsBase, Statistics, .Threads, JLD2, DataFrames

## Create data frame of means and stds for each summary statistic
μσ_df = Dict(S .=> [DataFrame(
            hcat([[L,t,
                mean(filter(isfinite,Data[(Data.PoreSize .== L) .& (Data.Day .== t),Sᵢ])),
                std(filter(isfinite,Data[(Data.PoreSize .== L) .& (Data.Day .== t),Sᵢ]))
                  ] for t ∈ T, L ∈ P]...)',[:PoreSize,:Day,:μ,:σ]) for Sᵢ ∈ S])

## Fit quadratic for each summary statistic
σ  = Dict()

# Pre-specify intercept as 10% of the maximum
σ_intercepts = Dict(S .=> [0.1*maximum(filter(isfinite,μσ_df[Sᵢ].σ)) for Sᵢ ∈ S])

for (i,Sᵢ) ∈ enumerate(S)

    # Generic form of quadratic
    f = (x,θ) -> θ[1] * x^2 + θ[2] * x + σ_intercepts[Sᵢ]
    θ = FitModel(f,[0.0,0.0],μσ_df[Sᵢ].μ,μσ_df[Sᵢ].σ)
    push!(σ, Sᵢ => μ -> max.(σ_intercepts[Sᵢ],f(μ,θ)))

end

## Save results
@save "Results/StandardDeviations.jld2" μσ_df σ σ_intercepts
