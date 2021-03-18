#=
#
#   LogLikelihood.jl
#
#   Computes log likelihood for the pore model given data and options.
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au
#       https://alexbrowning.me
#
=#

function LogLikelihood(θ,
                Data;
                β   = ones(0),
                i_β = ones(Int,0),
                Cᵢ  = [:S1_Coverage],
                σ   = μ -> 1.0*length(μ),
                L   = 300.0,
                t₀  = 1.0,
                kwargs...
            )

    # Full parameter set
    n_θ = length(θ) + length(β)
    i_θ = setdiff(1:n_θ,i_β)
    Θ   = zeros(n_θ)
    Θ[i_β] = β
    Θ[i_θ] = θ
    D,λ,K,α,τ,u₀ = Θ

    # Times for inference
    T₁  = T[T .> t₀]

    # Subset data
    Data₁ = Data[(Data.PoreSize .== L),:]

    # Solve PDE
    sol = SolvePDE(D,λ,K,α,T₁,L=L,u₀=u₀,t₀=t₀,kwargs...)

    # Initialise log likelihood
    LL = 0.0

    # Loop through time points
    for (i_T,t) ∈ enumerate(T₁)

        # Model prediction (mean)
        μₘ = SummaryStatistics(sol(t),L=L,τₐ=τ*K)

        # Loop through summary statistics
        for Sᵢ ∈ Cᵢ

            # Standard deviation
            σₘ = σ[Sᵢ](μₘ[Sᵢ])

            # Loop through observations, yₒ ∼ N(μₘ,σₘ)
            Yₒ = Data₁[Data₁.Day .== t, Sᵢ]
            for yᵢ ∈ filter(isfinite,Yₒ)

                LL += -log(σₘ * √(2π)) - 0.5 * (yᵢ - μₘ[Sᵢ])^2 / σₘ^2

            end

        end

    end

    return LL

end

function GetΘ(θ,β,i_β)

    n_θ = length(θ) + length(β)
    i_θ = setdiff(1:n_θ,i_β)
    Θ   = zeros(n_θ)
    Θ[i_β] = β
    Θ[i_θ] = θ

    return Θ

end
#
# function LogLikelihood(θ,
#                 Data;
#                 β   = ones(0),
#                 i_β = ones(Int,0),
#                 Cᵢ  = [:S1_Coverage],
#                 Σ   = diagm(ones(length(Cᵢ))),
#                 L   = 300.0,
#                 t₀  = 1.0,
#                 kwargs...
#             )
#
#     detΣ = det(Σ)
#     invΣ = inv(Σ)
#
#     # Full parameter set
#     n_θ = length(θ) + length(β)
#     i_θ = setdiff(1:n_θ,i_β)
#     Θ   = zeros(n_θ)
#     Θ[i_β] = β
#     Θ[i_θ] = θ
#     D,λ,K,α,τ,u₀ = Θ
#
#     # Times for inference
#     T₁  = T[T .> t₀]
#
#     # Subset data
#     Data₁ = Data[(Data.PoreSize .== L),:]
#
#     # Solve PDE
#     sol = SolvePDE(D,λ,K,α,T₁,L=L,u₀=u₀,t₀=t₀,kwargs...)
#
#     # Initialise log likelihood
#     LL = 0.0
#
#     # Loop through time points
#     for (i_T,t) ∈ enumerate(T₁)
#
#         # Model prediction (mean)
#         μₘ_all = SummaryStatistics(sol(t),L=L,τₐ=τ*K)
#         μₘ = [μₘ_all[Sᵢ] for Sᵢ ∈ Cᵢ]
#
#         Xₒ = hcat([Data₁[Data₁.Day .== t, Sᵢ] for Sᵢ ∈ Cᵢ]...)
#
#         # Loop through observations
#         for i = 1:size(Xₒ,1)
#
#             xₒ = Xₒ[i,:]
#             if all(isfinite.(xₒ))
#
#                 LL += -2/2 * log(2π) - 0.5 * log(detΣ) - 0.5 * (xₒ - μₘ)' * invΣ * (xₒ - μₘ)
#
#             end
#
#         end
#
#     end
#
#     return LL
#
# end
