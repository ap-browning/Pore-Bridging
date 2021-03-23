#=
#
#   Profiles2D.jl
#
#   Compute bivadiate profile likelihood functions.
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

# Summary statistic combinations to use
C = [:S1_Density, :S2_Coverage]

# Bounds to profile
#        D        λ       K       u₀
lb_p  = [10.0,    0.1,    0.0025,  1e-5]
ub_p  = [1500.0,  2.0,    0.0045,  0.002]

# Initialise storage
Profiles2D = [Dict() for i = 1:length(P)]

# Compute and time main result
t1 = @elapsed begin

    # Loop over C and P
    @threads for i_P = 1:4
        L = P[i_P]

        # Maximum likelihood estimate
        θ̂     = MLE_Individual[1][i_P][:θ̂]
        L̂     = MLE_Individual[1][i_P][:L̂]

        # Values to profile
        Ψ     = collect.(range.(lb,ub,length=profile_res_biv))

        ## Profile over 1 and 2
        bv_p̂₁ = zeros(length(Ψ[1]),length(Ψ[2]))

            # Loop over parameter 1 (start close to MLE)
            idx = findmin(abs.(Ψ[1] .- θ̂[1]))[2]

            # Initial guess for maximum with parameter 1 fixed is the MLE
            λ̂ = θ̂[[2,3,4]]

            # Above MLE
            for i = idx:length(Ψ[1])

                # Fix parameter 1
                β₂      = [β;Ψ[1][i]]
                i_β₂    = [i_β;1]

                # New likelihood function
                l    = θ -> LogLikelihood(θ,Data₁;β=β₂,i_β=i_β₂,L=L,Cᵢ=C,σ=σ,t₀=t₀)

                # Update guess for maximum with parameter 1 fixed
                λ̂,   = Maximise(l,λ̂,:LN_BOBYQA,lb=lb[[2,3,4]],ub=ub[[2,3,4]],ftol_abs=ftol_abs)

                # Profile
                bv_p̂₁[i,:] = Profile(l,Ψ[2],1,θ̂=λ̂,lb=lb[[3,4]],ub=ub[[3,4]],ftol_abs=ftol_abs)

            end

            # Below MLE
            λ̂ = θ̂[[2,3,4]]
            for i = idx-1:-1:1

                # Fix parameter 1
                β₂      = [β;Ψ[1][i]]
                i_β₂    = [i_β;1]

                # New likelihood function
                l    = θ -> LogLikelihood(θ,Data₁;β=β₂,i_β=i_β₂,L=L,Cᵢ=C,σ=σ,t₀=t₀)

                # Update guess for maximum with parameter 1 fixed
                λ̂,   = Maximise(l,λ̂,:LN_BOBYQA,lb=lb[[2,3,4]],ub=ub[[2,3,4]],ftol_abs=ftol_abs)

                # Profile
                bv_p̂₁[i,:] = Profile(l,Ψ[2],1,θ̂=λ̂,lb=lb[[3,4]],ub=ub[[3,4]],ftol_abs=ftol_abs)

            end

        ## Profile over 2 and 3
        bv_p̂₂ = zeros(length(Ψ[2]),length(Ψ[3]))

            # Loop over parameter 2 (start close to MLE)
            idx = findmin(abs.(Ψ[2] .- θ̂[2]))[2]

            # Initial guess for maximum with parameter 1 fixed is the MLE
            λ̂ = θ̂[[1,3,4]]

            # Above MLE
            for i = idx:length(Ψ[2])

                # Fix parameter 2
                β₂      = [β;Ψ[2][i]]
                i_β₂    = [i_β;2]

                # New likelihood function
                l    = θ -> LogLikelihood(θ,Data₁;β=β₂,i_β=i_β₂,L=L,Cᵢ=C,σ=σ,t₀=t₀)

                # Update guess for maximum with parameter 2 fixed
                λ̂,   = Maximise(l,λ̂,:LN_BOBYQA,lb=lb[[1,3,4]],ub=ub[[1,3,4]],ftol_abs=ftol_abs)

                # Profile
                bv_p̂₂[i,:] = Profile(l,Ψ[3],2,θ̂=λ̂,lb=lb[[1,4]],ub=ub[[1,4]],ftol_abs=ftol_abs)

            end

            # Below MLE
            λ̂ = θ̂[[1,3,4]]
            for i = idx-1:-1:1

                # Fix parameter 2
                β₂      = [β;Ψ[2][i]]
                i_β₂    = [i_β;2]

                # New likelihood function
                l    = θ -> LogLikelihood(θ,Data₁;β=β₂,i_β=i_β₂,L=L,Cᵢ=C,σ=σ,t₀=t₀)

                # Update guess for maximum with parameter 2 fixed
                λ̂,   = Maximise(l,λ̂,:LN_BOBYQA,lb=lb[[1,3,4]],ub=ub[[1,3,4]],ftol_abs=ftol_abs)

                # Profile
                bv_p̂₂[i,:] = Profile(l,Ψ[3],2,θ̂=λ̂,lb=lb[[1,4]],ub=ub[[1,4]],ftol_abs=ftol_abs)

            end

        # Save results
        push!(Profiles2D[i_P], :Ψ  => Ψ)
        push!(Profiles2D[i_P], :bv_p̂₁ => bv_p̂₁)
        push!(Profiles2D[i_P], :bv_p̂₂ => bv_p̂₂)

    end

end

# Save settings
Profiles2D_Settings = Dict(
    :β  => β,  :i_β => i_β, :lb  => lb, :ub => ub,
    :θ₀ => θ₀, :t₀  => t₀,  :C   => C,
    :ftol_abs => ftol_abs,  :walltime => t1,
    :profile_res_biv => profile_res_biv,
    :MLE_results => MLE_Individual, :MLE_settings => MLE_Individual_Settings)

# Save
@save "Results/Saved/Profiles2D.jld2" Profiles2D Profiles2D_Settings
