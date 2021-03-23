
using Pores
using JLD2, PyPlot, StatsBase
include("../Preferences.jl")

# Create figure
figS4a,axs = subplots(1,4,figsize=(6.7,1.6))

# Load results
@load "Results/Supplementary/Saved/All400_MLE_Individual.jld2"

# Options
LineStyles = ["-","--",":"]
Tfine      = collect(range(4.0,28.0,length=200))

# 400μm pores
i_P,L = 2,400.0

# Loop over S
for (i_S,Sᵢ) ∈ enumerate(S)

    # Current axis
    ax    = axs[i_S]

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
    for (i_C,Cᵢ) ∈ enumerate(All400_MLE_Individual_Settings[:C])

        # MLE
        D,λ,K,α,τ,u₀ = GetΘ(All400_MLE_Individual[i_C][i_P][:θ̂],All400_MLE_Individual_Settings[:β],All400_MLE_Individual_Settings[:i_β])

        # Solve and summarise PDE
        sol = SolvePDE(D,λ,K,α,Tfine,L=L,t₀=All400_MLE_Individual_Settings[:t₀],u₀=u₀)
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
    ax.set_ylabel(Snms[Sᵢ])
    ax.set_xlabel("Time (d)")

end

NumberPlots!(axs)
plt.tight_layout(pad=0.0,h_pad=1.0,w_pad=1.0)
savefig("Figures/Supplementary/Saved/Fig-S4a-All400_Fits.pdf")
display(figS4a)
