#=
#
#   Profile.jl
#
#   Perform bounded maximisation using methods in NLopt.jl
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au
#       https://alexbrowning.me
#
=#
@doc """
    Profile()

Maximises function `fun` (of the single vector variable `θ0`) within bounds
    specified by `lb` (lower) and `ub` (upper). Uses `NLopt.jl` routine given by
    `method` where the tolerence is given by `f_tol`.

Inputs:\n
    `fun`       - Function to maximise, of variable `θ`, fun(θ) : R^N → R
    `θ0`        - Initial guess
    `lb`        - lower bounds (default: zeros)
    `ub`        - upper bounds (default: Inf)
    `ftol_rel`  - tolerance (default: 0.0)
    `ftol_abs`  - tolerance (default: 0.0)
    `method`    - method to use (default: COBYLA)

Outputs:\n
    `θ̂` - argmax(fun(θ))
    `L̂` - max(fun(θ))
"""
# CLEVER PROFILE
function Profile(fun::Function,ψ::Array{Float64,1},i_ψ::Int64,
                 method=:LN_BOBYQA;
                 θ̂::Array{Float64,1},
                 kwargs...)

    # Initialise output
    p̂ = zeros(length(ψ))

    # Indices of nuisance parameters
    n_θ = length(θ̂)
    i_λ = setdiff(1:n_θ,i_ψ)

    # Parameter map function
    function LL_λ(ψ,λ)

        θ = zeros(n_θ)
        θ[i_ψ] = ψ
        θ[i_λ] = λ

        return fun(θ)

    end

    # Where the maximum should be (start near the MLE)
    ψ̂   = θ̂[i_ψ]
    above_idx = ψ[end] > ψ̂ ? findfirst(ψ .> ψ̂) : length(ψ) + 1

    # Above the MLE
    λ₀  = θ̂[i_λ]
    for i = above_idx:length(ψ)

        λ₀,p̂[i] = Maximise(λ -> LL_λ(ψ[i],λ),λ₀,method;kwargs...)

    end

    # Below the MLE
    λ₀  = θ̂[i_λ]
    for i = above_idx-1:-1:1

        λ₀,p̂[i] = Maximise(λ -> LL_λ(ψ[i],λ),λ₀,method;kwargs...)

    end

    return p̂

end

function ConfidenceInterval(Ψ,p̂)

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
