
using Pores
using JLD2, PyPlot, StatsBase
include("../Preferences.jl")

# Create figure
figS8,axs = subplots(2,4,figsize=(6.7,3))

# Load results
@load "Results/Supplementary/Saved/Fisher_MLE_Individual.jld2"

# Options
LineStyles = ["-","--",":"]
Tfine      = collect(range(4.0,28.0,length=200))

# Correlations to plot
ρ = [[:S2_Coverage,:S1_Density], [:S2_Coverage, :S3_Circularity]]

# Loop over P and correlations
for (i_P,L) ∈ enumerate(P), (i_ρ,ρᵢ) ∈ enumerate(ρ)

    # Current axis
    ax = axs[i_ρ,i_P]

    # Plot data
    ax.plot(Data[Data.PoreSize .== L,ρᵢ[1]],Data[Data.PoreSize .== L,ρᵢ[2]],PoreSizeSymb[i_P],ms=PoreSizeMkSz[i_P],alpha=0.15,c=PoreSizeColsLight[i_P])

    # Plot solution where all pore sizes considered simultaneously
    
        # MLE
        D,λ,K,α,τ,u₀ = GetΘ(Fisher_MLE_Combined[:θ̂],Fisher_MLE_Combined_Settings[:β],Fisher_MLE_Combined_Settings[:i_β])[[1,2,3,4,5,5+i_P]]
    
        # Solve and summarise PDE
        sol = SolvePDE(D,λ,K,α,Tfine,L=L,t₀=Fisher_MLE_Individual_Settings[:t₀],u₀=u₀)
        Xplot = zeros(size(Tfine))
        Yplot = zeros(size(Tfine))
        for (i,t) ∈ enumerate(Tfine)
            Xplot[i] = SummaryStatistics(sol(t),L=L,τₐ=τ*K)[ρᵢ[1]]
            Yplot[i] = SummaryStatistics(sol(t),L=L,τₐ=τ*K)[ρᵢ[2]]
        end
        ax.plot(Xplot,Yplot,"-",color="#333333",lw=1)

    # Plot solution for each set of summary statistics
    for (i_C,Cᵢ) ∈ enumerate(Fisher_MLE_Individual_Settings[:C])

        # MLE
        D,λ,K,α,τ,u₀ = GetΘ(Fisher_MLE_Individual[i_C][i_P][:θ̂],Fisher_MLE_Individual_Settings[:β],Fisher_MLE_Individual_Settings[:i_β])

        # Solve and summarise PDE
        sol = SolvePDE(D,λ,K,α,Tfine,L=L,t₀=Fisher_MLE_Individual_Settings[:t₀],u₀=u₀)
        Xplot = zeros(size(Tfine))
        Yplot = zeros(size(Tfine))
        for (i,t) ∈ enumerate(Tfine)
            Xplot[i] = SummaryStatistics(sol(t),L=L,τₐ=τ*K)[ρᵢ[1]]
            Yplot[i] = SummaryStatistics(sol(t),L=L,τₐ=τ*K)[ρᵢ[2]]
        end
        ax.plot(Xplot,Yplot,LineStyles[i_C],color=PoreSizeCols[i_P],lw=1)

    end

    # Plot "circle before close" marker
    if i_ρ == 2
        ax.plot(1.0 - π/4,1.0,"k.")
    end

    # Style axes
    ax.set_xlim(Slim[ρᵢ[1]])
    ax.set_ylim(Slim[ρᵢ[2]])
    ax.set_xticks(Stks[ρᵢ[1]])
    ax.set_yticks(Stks[ρᵢ[2]])
    if i_P == 1
        ax.set_ylabel(Snms[ρᵢ[2]])
    else
        ax.set_yticklabels([])
    end
    if i_ρ == 1
        IntL = Int(L)
        ax.set_title("$IntL μm")
    end
    if i_ρ == length(ρ)
        ax.set_xlabel(Snms[ρᵢ[1]])
    else
        ax.set_xticklabels([])
    end

end

NumberPlots!(axs)
plt.tight_layout(pad=0.0,h_pad=1.0,w_pad=1.0)
savefig("Figures/Supplementary/Saved/Fig-S8-Fisher_Correlations.pdf")
display(figS8)
