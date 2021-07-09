
using Pores
using JLD2, PyPlot, StatsBase
include("../Preferences.jl")

# Create figure
figS12,axs = subplots(4,4,figsize=(6.7,5))

# Load results
@load "Results/Saved/MLE_Combined.jld2"

# Second grid spacing to try
N₂ = 151

# Options
Tfine      = collect(range(4.0,28.0,length=200))

# Loop over S and P
for (i_P,L) ∈ enumerate(P)

    # MLE
    D,λ,K,α,τ,u₀ = GetΘ(MLE_Combined[:θ̂],MLE_Combined_Settings[:β],MLE_Combined_Settings[:i_β])[[1,2,3,4,5,5+i_P]]

    # Solve PDE (grid 1)
    sol1 = SolvePDE(D,λ,K,α,Tfine,L=L,t₀=MLE_Combined_Settings[:t₀],u₀=u₀)

    # Solve PDE (grid 2)
    sol2 = SolvePDE(D,λ,K,α,Tfine,L=L,t₀=MLE_Combined_Settings[:t₀],u₀=u₀,N=N₂)

    for (i_S,Sᵢ) ∈ enumerate(S)
        
        # Current axis
        ax    = axs[i_S,i_P]

        # Summarise PDE (grid 1)
        Xplot = zeros(size(Tfine))
        for (i,t) ∈ enumerate(Tfine)
            Xplot[i] = SummaryStatistics(sol1(t),L=L,τₐ=τ*K)[Sᵢ]
        end
        ax.plot(Tfine,Xplot,"-",color=PoreSizeCols[i_P],lw=1)

        # Summarise PDE (grid 2)
        Xplot = zeros(size(Tfine))
        for (i,t) ∈ enumerate(Tfine)
            Xplot[i] = SummaryStatistics(sol2(t),L=L,τₐ=τ*K)[Sᵢ]
        end
        ax.plot(Tfine,Xplot,":k",lw=1)

        # Style axes
        ax.set_xlim([2.0,30.0])
        ax.set_xticks(T)
        ax.set_ylim(Slim[Sᵢ])
        ax.set_yticks(Stks[Sᵢ])
        if i_P == 1
            ax.set_ylabel(Snms[Sᵢ])
        else
            ax.set_yticklabels([])
        end
        if i_S == 1
            IntL = Int(L)
            ax.set_title("$IntL μm")
        end
        if i_S == length(S)
            ax.set_xlabel("Time (d)")
        else
            ax.set_xticklabels([])
        end
    end

end

NumberPlots!(axs)
plt.tight_layout(pad=0.0,h_pad=1.0,w_pad=1.0)
savefig("Figures/Supplementary/Saved/Fig-S12-GridDependence2.pdf")
display(figS12)