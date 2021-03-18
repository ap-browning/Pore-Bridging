#=

    MLE.jl

    Compute maximum likelihood estimates for the Porous-Fisher model with
    quadratic noise function. Repeat analysis using the following summary
    statistic combinations.

        1. S1_Density and S2_Coverage
        2. Only S1_Density
        3. S2_Coverage (and end-time S1_Density)

=#

# Load module
using PoresIdentifiability
using StatsBase, Statistics, .Threads, JLD2, TimerOutputs

# Load options
include("Options.jl")

# Load noise model
include("StandardDeviations.jl")

## INDIVIDUAL (PER PORE SIZE)

# Summary statistic combinations to use
C = [[:S1_Density, :S2_Coverage], [:S1_Density], [:S2_Coverage]]

# Initialise storage
MLE_Individual = [[Dict() for i = 1:length(P)] for i = 1:length(C)]

# Loop over C and P
t1 = @elapsed begin

    @threads for i_P = 1:4
        L = P[i_P]
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
            θ̂,L̂   = Maximise(l,θ₀,:GN_DIRECT,lb=lb,ub=ub,maxtime=maxtime)
            θ̂,L̂   = Maximise(l,θ̂,:LN_BOBYQA,lb=lb,ub=ub,ftol_abs=ftol_abs)

            # Save results
            push!(MLE_Individual[i_C][i_P], :θ̂  => θ̂)
            push!(MLE_Individual[i_C][i_P], :L̂  => L̂)

        end

    end

end

# Save settings
MLE_Individual_Settings = Dict(
    :β  => β,  :i_β => i_β, :lb  => lb, :ub => ub,
    :θ₀ => θ₀, :t₀  => t₀,  :C   => C,
    :ftol_abs => ftol_abs,  :maxtime => maxtime,
    :walltime => t1)

# Save
@save "Results/MLE_Individual.jld2" MLE_Individual MLE_Individual_Settings

## COMBINED (ALL PORE SIZES)

# Summary statistic combinations to use
C   = [:S1_Density,:S2_Coverage]

# Adjust parameter space to allow initial condition to vary
    # Parameter bounds
    #                |---------- u₀ ----------|
    #                 300    400    500    600
    lb₁  = [lb[1:3]; [lb[4], lb[4], lb[4], lb[4]] ]
    ub₁  = [ub[1:3]; [ub[4], ub[4], ub[4], ub[4]] ]

    # Initial guess (algorithm actually independent of this choice)
    θ₀  = lb₁ + rand(length(lb₁)) .* (ub₁ - lb₁)

# Initialise storage
MLE_Combined = Dict()

# Calculate combined MLE
t2 = @elapsed begin

    # Construct log likelihood function
    l = θ -> LogLikelihood(θ[[1,2,3,4]],Data₁;β=β,i_β=i_β,L=300.0,Cᵢ=C,σ=σ,t₀=t₀) +
             LogLikelihood(θ[[1,2,3,5]],Data₁;β=β,i_β=i_β,L=400.0,Cᵢ=C,σ=σ,t₀=t₀) +
             LogLikelihood(θ[[1,2,3,6]],Data₁;β=β,i_β=i_β,L=500.0,Cᵢ=C,σ=σ,t₀=t₀) +
             LogLikelihood(θ[[1,2,3,7]],Data₁;β=β,i_β=i_β,L=600.0,Cᵢ=C,σ=σ,t₀=t₀)

    # Find maximum likelihood estimate
    θ̂,L̂   = Maximise(l,θ₀,:GN_DIRECT,lb=lb₁,ub=ub₁,maxtime=maxtime)
    θ̂,L̂   = Maximise(l,θ̂,:LN_BOBYQA,lb=lb₁,ub=ub₁,ftol_abs=ftol_abs)

    # Save results
    push!(MLE_Combined, :θ̂  => θ̂)
    push!(MLE_Combined, :L̂  => L̂)

end

# Save settings
MLE_Combined_Settings = Dict(
    :β  => β,  :i_β => i_β, :lb₁  => lb₁, :ub₁ => ub₁,
    :θ₀ => θ₀, :t₀  => t₀,  :C   => C,
    :ftol_abs => ftol_abs,  :maxtime => maxtime,
    :walltime => t2)

# Save
@save "Results/MLE_Combined.jld2" MLE_Combined MLE_Combined_Settings
