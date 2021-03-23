#=
#
#   Fisher_Profiles.jl
#
#   Repeat results in Profiles.jl for Fisher's equation (α = 0.0)
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

# Load modules
using Pores
using StatsBase, Statistics, .Threads, JLD2, TimerOutputs

# Load options
include("Fisher_Options.jl")

# Load noise model
include("../StandardDeviations.jl")

# Load maximum likelihood estimates
@load "Results/Supplementary/Saved/Fisher_MLE_Individual.jld2"

# Summary statistic combinations to use
C = [[:S1_Density, :S2_Coverage]]

# Initialise storage
Fisher_Profiles_Individual = [[Dict() for i = 1:length(P)] for i = 1:length(C)]

# Compute and time main result
t1 = @elapsed begin

    # Loop over C and P
    @time for (i_C,Cᵢ) ∈ enumerate(C), (i_P,L) ∈ enumerate(P)

        # Construct log likelihood function
        l     = θ -> LogLikelihood(θ,Data₁;β=β,i_β=i_β,L=L,Cᵢ=Cᵢ,σ=σ,t₀=t₀)

        # Values to profile
        Ψ     = collect.(range.(lb,ub,length=profile_res_uni))

        # Maximum likelihood estimate
        θ̂     = Fisher_MLE_Individual[i_C][i_P][:θ̂]
        L̂     = Fisher_MLE_Individual[i_C][i_P][:L̂]

        # Profile over all parameters
        p̂    = Array{Any,1}(undef,4)
        @threads for i_ψ = 1:4
            i_λ     = setdiff(1:length(θ̂),i_ψ)
            p̂[i_ψ]  = Profile(l,Ψ[i_ψ],i_ψ,θ̂=θ̂,lb=lb[i_λ],ub=ub[i_λ],ftol_abs=ftol_abs)
        end

        # Save results
        push!(Fisher_Profiles_Individual[i_C][i_P], :θ̂  => θ̂)
        push!(Fisher_Profiles_Individual[i_C][i_P], :L̂  => L̂)
        push!(Fisher_Profiles_Individual[i_C][i_P], :p̂  => p̂)
        push!(Fisher_Profiles_Individual[i_C][i_P], :Ψ  => Ψ)

    end

end

# Save settings
Fisher_Profiles_Individual_Settings = Dict(
    :β  => β,  :i_β => i_β, :lb  => lb, :ub => ub,
    :θ₀ => θ₀, :t₀  => t₀,  :C   => C,
    :ftol_abs => ftol_abs,  :walltime => t1,
    :profile_res_uni => profile_res_uni,
    :MLE_results => Fisher_MLE_Individual, :MLE_settings => Fisher_MLE_Individual_Settings)

# Save
@save "Results/Supplementary/Saved/Fisher_Profiles_Individual.jld2" Fisher_Profiles_Individual Fisher_Profiles_Individual_Settings
