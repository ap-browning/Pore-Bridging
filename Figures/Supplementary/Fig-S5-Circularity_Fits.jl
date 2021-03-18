
using PoresIdentifiability
using JLD2, PyPlot, StatsBase
include("../Preferences.jl")

# Create figure
figS5,axs = subplots(4,4,figsize=(6.7,5))

# Load results
@load "Results/Supplementary/Circularity_MLE_Individual.jld2"

# Options
LineStyles = ["-."]
Tfine      = collect(range(4.0,28.0,length=200))

# Loop over S and P
for (i_S,Sᵢ) ∈ enumerate(S), (i_P,L) ∈ enumerate(P)

    # Current axis
    ax    = axs[i_S,i_P]

    # Plot data
    Dplot = [filter(isfinite,Data[(Data.PoreSize .== L) .& (Data.Day .== t), Sᵢ]) for t ∈ T]
    Tplot = T[(!).(isempty.(Dplot))]
    Dplot = Dplot[(!).(isempty.(Dplot))]
    vp    = ax.violinplot(Dplot,Tplot,widths=2)
    for pc = vp["bodies"]
        pc.set_facecolor(PoreSizeCols[i_P])
        pc.linewidths = 0
    end
    vp["cbars"].set_linewidth(0)
    vp["cmins"].set_linewidth(0)
    vp["cmaxes"].set_linewidth(0)

    # Plot solution for each set of summary statistics
    for (i_C,Cᵢ) ∈ enumerate(Circularity_MLE_Individual_Settings[:C])

        # MLE
        D,λ,K,α,τ,u₀ = GetΘ(Circularity_MLE_Individual[i_C][i_P][:θ̂],Circularity_MLE_Individual_Settings[:β],Circularity_MLE_Individual_Settings[:i_β])

        # Solve and summarise PDE
        sol = SolvePDE(D,λ,K,α,Tfine,L=L,t₀=Circularity_MLE_Individual_Settings[:t₀],u₀=u₀)
        Xplot = zeros(size(Tfine))
        for (i,t) ∈ enumerate(Tfine)
            Xplot[i] = SummaryStatistics(sol(t),L=L,τₐ=τ*K)[Sᵢ]
        end
        ax.plot(Tfine,Xplot,LineStyles[i_C],color=PoreSizeCols[i_P],lw=1)

    end

    # Style axes
    ax.set_xlim([4-2.0,28+2.0])
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
        ax.set_xlabel("Day")
    else
        ax.set_xticklabels([])
    end

end

NumberPlots!(axs)
plt.tight_layout(pad=0.0,h_pad=1.0,w_pad=1.0)
savefig("Figures/Supplementary/Fig-S5-Circularity_Fits.pdf")
display(figS5)
