#=
#
#   LogLikelihood.jl
#
#   Computes log likelihood for the pore model given data and options.
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

function LogLikelihood(θ,                   # Parameters
                Data;                       # Dataframe used for inference
                β   = ones(0),              # Fixed parameters
                i_β = ones(Int,0),          # Indices of fixed parameters
                Cᵢ  = [:S1_Coverage],       # Summary statistics included in likelihood function
                σ   = μ -> 1.0*length(μ),   # Standard deviation as function of mean (should be per summary statistic)
                L   = 300.0,                # Domain size
                t₀  = 1.0,                  # Time at which to apply the initial condition
                kwargs...                   # Extra arguments to pass to PDE solver
            )

    # Full parameter set, Θ = (θ,β) (need to reorder)
    Θ = GetΘ(θ,β,i_β)
    D,λ,K,α,τ,u₀ = Θ

    # Times included in inference
    T₁  = T[T .> t₀]

    # Subset data
    Data₁ = Data[(Data.PoreSize .== L),:]

    # Solve PDE
    sol = SolvePDE(D,λ,K,α,T₁,L=L,u₀=u₀,t₀=t₀,kwargs...)

    # Initialise log likelihood
    LL = 0.0

    # Loop through time points
    for (i_T,t) ∈ enumerate(T₁)

        # Summarise PDE solution at time t as thge model prediction (mean)
        μₘ = SummaryStatistics(sol(t),L=L,τₐ=τ*K)

        # Loop through summary statistics that are included
        for Sᵢ ∈ Cᵢ

            # Standard deviation for each summary statistic
            σₘ = σ[Sᵢ](μₘ[Sᵢ])

            # Loop through observations, yₒ ∼ N(μₘ,σₘ) and add to log-likelihood
            Yₒ = Data₁[Data₁.Day .== t, Sᵢ]
            for yᵢ ∈ filter(isfinite,Yₒ)

                LL += -log(σₘ * √(2π)) - 0.5 * (yᵢ - μₘ[Sᵢ])^2 / σₘ^2

            end # yᵢ ∈ filter(isfinite,Yₒ)

        end # for Sᵢ ∈ Cᵢ

    end # (i_T,t) ∈ enumerate(T₁)

    # Return log-likelihood function
    return LL

end

# Sub-routine that converts partitioned parameter space (θ,β) into full parameter space Θ
function GetΘ(θ,β,i_β)

    n_θ = length(θ) + length(β)
    i_θ = setdiff(1:n_θ,i_β)
    Θ   = zeros(n_θ)
    Θ[i_β] = β
    Θ[i_θ] = θ

    return Θ

end