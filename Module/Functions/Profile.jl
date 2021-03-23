#=
#
#   Profile.jl
#
#   Compute profile likelihood in the neighbourhood of the MLE θ̂
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

# Example usage:
#   Suppose we have θ = (D,λ,K) and we wish to profile λ between 0.0 and 2.0
#
#       fun = θ -> LogLikelihood(θ,...)
#       ψ   = collect(range(0.0,2.0,length=50))
#       i_ψ = 2
#       θ̂,  = Maximise(fun,...)
#       
#       p̂   = Profile(fun,ψ,i_ψ,θ̂=θ̂)

function Profile(fun::Function,
                 ψ::Array{Float64,1},   # Vector of values of profiled parameter
                 i_ψ::Int64,            # Index of parameter to profile
                 method=:LN_BOBYQA;
                 θ̂::Array{Float64,1},
                 kwargs...)

    # Initialise output
    p̂ = zeros(length(ψ))

    # Indices of nuisance parameters, here θ = (ψ,λ)
    n_θ = length(θ̂)
    i_λ = setdiff(1:n_θ,i_ψ)

    # Call fun(θ) using ψ and λ
    function LL_λ(ψ,λ)

        θ = zeros(n_θ)
        θ[i_ψ] = ψ
        θ[i_λ] = λ

        return fun(θ)

    end

    # Start close to the MLE
    ψ̂   = θ̂[i_ψ]
    above_idx = ψ[end] > ψ̂ ? findfirst(ψ .> ψ̂) : length(ψ) + 1

    # Profile above the MLE
    λ₀  = θ̂[i_λ]
    for i = above_idx:length(ψ)

        # Optimise for ψ fixed
        λ₀,p̂[i] = Maximise(λ -> LL_λ(ψ[i],λ),λ₀,method;kwargs...)

    end

    # Profile below the MLE
    λ₀  = θ̂[i_λ]
    for i = above_idx-1:-1:1

        # Optimise for ψ fixed
        λ₀,p̂[i] = Maximise(λ -> LL_λ(ψ[i],λ),λ₀,method;kwargs...)

    end

    # Return profile likelihood
    return p̂

end

# Compute approximate confidence interval given profile likelihood (only works for "nice" functions)
function ConfidenceInterval(Ψ,p̂)

    # Normalise, here max(p) = 0.0
    p = p̂ .- maximum(p̂)

    # Calculate 95% CI
    lwr,upr = Ψ[[1,end]]
    if p[1] < -1.92
        idx = findfirst(p .> -1.92) - 1
        y₁,y₂ = p[[idx,idx+1]]
        x₁,x₂ = Ψ[[idx,idx+1]]
        lwr   = x₁ + (-1.92 - y₁) * (x₂ - x₁) / (y₂ - y₁)
        display(lwr)
    end
    if p[end] < -1.92
        idx = findlast(p .> -1.92)
        y₁,y₂ = p[[idx,idx+1]]
        x₁,x₂ = Ψ[[idx,idx+1]]
        upr   = x₁ + (-1.92 - y₁) * (x₂ - x₁) / (y₂ - y₁)
    end

    return lwr,upr

end
