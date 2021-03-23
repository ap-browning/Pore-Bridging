#=
#
#   Profiles.jl
#
#   Compute univariate profile likelihood functions for each unknown parameter.
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
include("Options.jl")

# Load noise model
include("StandardDeviations.jl")

# Load maximum likelihood estimates
@load "Results/Saved/MLE_Individual.jld2"
@load "Results/Saved/MLE_Combined.jld2"

################################################
## INDIVIDUAL (PER PORE SIZE)

# Summary statistic combinations to use
C = [[:S1_Density, :S2_Coverage], [:S1_Density]]

# Initialise storage
Profiles_Individual = [[Dict() for i = 1:length(P)] for i = 1:length(C)]

# Compute and time main result
t1 = @elapsed begin

    # Loop over C and P
    @time for (i_C,Cᵢ) ∈ enumerate(C), (i_P,L) ∈ enumerate(P)

        # Construct log likelihood function
        l     = θ -> LogLikelihood(θ,Data₁;β=β,i_β=i_β,L=L,Cᵢ=Cᵢ,σ=σ,t₀=t₀)

        # Values to profile
        Ψ     = collect.(range.(lb,ub,length=profile_res_uni))

        # Maximum likelihood estimate
        θ̂     = MLE_Individual[i_C][i_P][:θ̂]
        L̂     = MLE_Individual[i_C][i_P][:L̂]

        # Profile over all parameters
        p̂    = Array{Any,1}(undef,4)
        @threads for i_ψ = 1:4
            i_λ     = setdiff(1:length(θ̂),i_ψ)
            p̂[i_ψ]  = Profile(l,Ψ[i_ψ],i_ψ,θ̂=θ̂,lb=lb[i_λ],ub=ub[i_λ],ftol_abs=ftol_abs)
        end

        # Save results
        push!(Profiles_Individual[i_C][i_P], :θ̂  => θ̂)
        push!(Profiles_Individual[i_C][i_P], :L̂  => L̂)
        push!(Profiles_Individual[i_C][i_P], :p̂  => p̂)
        push!(Profiles_Individual[i_C][i_P], :Ψ  => Ψ)

    end

end

# Save settings
Profiles_Individual_Settings = Dict(
    :β  => β,  :i_β => i_β, :lb  => lb, :ub => ub,
    :θ₀ => θ₀, :t₀  => t₀,  :C   => C,
    :ftol_abs => ftol_abs,  :walltime => t1,
    :profile_res_uni => profile_res_uni,
    :MLE_results => MLE_Individual, :MLE_settings => MLE_Individual_Settings)

# Save
@save "Results/Saved/Profiles_Individual.jld2" Profiles_Individual Profiles_Individual_Settings

################################################
## COMBINED (ALL PORE SIZES)

# Summary statistic combinations to use
C   = [:S1_Density,:S2_Coverage]

# Adjust parameter space to allow initial condition to vary
    # Parameter bounds
    #                |---------- u₀ ----------|
    #                 300    400    500    600
    lb₁  = [lb[1:3]; [lb[4], lb[4], lb[4], lb[4]] ]
    ub₁  = [ub[1:3]; [ub[4], ub[4], ub[4], ub[4]] ]

# Initialise storage
Profiles_Combined = Dict()

# Calculate combined MLE
t2 = @elapsed begin

    # Construct log likelihood function
    l = θ -> LogLikelihood(θ[[1,2,3,4]],Data₁;β=β,i_β=i_β,L=300.0,Cᵢ=C,σ=σ,t₀=t₀) +
             LogLikelihood(θ[[1,2,3,5]],Data₁;β=β,i_β=i_β,L=400.0,Cᵢ=C,σ=σ,t₀=t₀) +
             LogLikelihood(θ[[1,2,3,6]],Data₁;β=β,i_β=i_β,L=500.0,Cᵢ=C,σ=σ,t₀=t₀) +
             LogLikelihood(θ[[1,2,3,7]],Data₁;β=β,i_β=i_β,L=600.0,Cᵢ=C,σ=σ,t₀=t₀)

    # Find maximum likelihood estimate
    θ̂ = MLE_Combined[:θ̂]
    L̂ = MLE_Combined[:L̂]

    # Values to profile
    Ψ     = collect.(range.(lb₁,ub₁,length=profile_res_uni))

    # Profile over all parameters
    p̂     = Array{Any,1}(undef,4)
    @threads for i_ψ = 1:3
        i_λ     = setdiff(1:length(θ̂),i_ψ)
        p̂[i_ψ]  = Profile(l,Ψ[i_ψ],i_ψ,θ̂=θ̂,lb=lb₁[i_λ],ub=ub₁[i_λ],ftol_abs=ftol_abs)
    end

    # Save results
    push!(Profiles_Combined, :θ̂  => θ̂)
    push!(Profiles_Combined, :L̂  => L̂)
    push!(Profiles_Combined, :p̂  => p̂)
    push!(Profiles_Combined, :Ψ  => Ψ)

end

# Save settings
Profiles_Combined_Settings = Dict(
    :β  => β,  :i_β => i_β, :lb₁  => lb₁, :ub₁ => ub₁,
    :θ₀ => θ₀, :t₀  => t₀,  :C   => C,
    :ftol_abs => ftol_abs,  :walltime => t2,
    :profile_res_uni => profile_res_uni,
    :MLE_results => MLE_Combined, :MLE_settings => MLE_Combined_Settings)

# Save
@save "Results/Saved/Profiles_Combined.jld2" Profiles_Combined Profiles_Combined_Settings
