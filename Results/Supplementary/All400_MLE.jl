#=
#
#   All400_MLE.jl
#
#   Repeat results in MLE.jl for the 400 μm pores, using all data.
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
include("All400_Options.jl")

# Load noise model
include("../StandardDeviations.jl")

# Summary statistic combinations to use
C = [[:S1_Density, :S2_Coverage], [:S1_Density], [:S2_Coverage]]

# Initialise storage
All400_MLE_Individual = [[Dict() for i = 1:length(P)] for i = 1:length(C)]

## Compute and time main result
t1 = @elapsed begin

    # Loop over P (skip all but 400μm pores)
    for i_P = 1:4

        if i_P != 2
            continue
        end
        L = P[i_P]

        # Loop over C
        for (i_C,Cᵢ) ∈ enumerate(C)

            # Construct log likelihood function
            l     = θ -> LogLikelihood(θ,Data₁;β=β,i_β=i_β,L=L,Cᵢ=Cᵢ,σ=σ,t₀=t₀)

            # Include late-time density data if looking at only coverage
            if Cᵢ == [:S2_Coverage]
                Data₂ = copy(Data₁)
                Data₂[Data₂.Day .< T[end],:S1_Density] .= NaN
                l     = θ -> LogLikelihood(θ,Data₂;β=β,i_β=i_β,L=L,Cᵢ=[:S1_Density,:S2_Coverage],σ=σ,t₀=t₀)
            end

            # Find maximum likelihood estimate
            θ̂,L̂   = Maximise(l,θ₀,:GN_DIRECT,lb=lb,ub=ub,maxtime=maxtime)       # Global routine
            θ̂,L̂   = Maximise(l,θ̂,:LN_BOBYQA,lb=lb,ub=ub,ftol_abs=ftol_abs)      # Local routine

            # Save results
            push!(All400_MLE_Individual[i_C][i_P], :θ̂  => θ̂)
            push!(All400_MLE_Individual[i_C][i_P], :L̂  => L̂)

        end

    end

end

# Save settings
All400_MLE_Individual_Settings = Dict(
    :β  => β,  :i_β => i_β, :lb  => lb, :ub => ub,
    :θ₀ => θ₀, :t₀  => t₀,  :C   => C,
    :ftol_abs => ftol_abs,  :maxtime => maxtime,
    :walltime => t1)

# Save
@save "Results/Supplementary/Saved/All400_MLE_Individual.jld2" All400_MLE_Individual All400_MLE_Individual_Settings
