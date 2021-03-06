
using Pores
using JLD2, PyPlot, StatsBase
include("Preferences.jl")

# Create figure
fig4,axs = subplots(4,4,figsize=(6.7,5))

# Load results
@load "Results/Saved/MLE_Individual.jld2"
@load "Results/Saved/MLE_Combined.jld2"

# Options
LineStyles = ["-","--",":"]
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

    # Plot solution where all pore sizes considered simultaneously

        #MLE
        D,λ,K,α,τ,u₀ = GetΘ(MLE_Combined[:θ̂],MLE_Combined_Settings[:β],MLE_Combined_Settings[:i_β])[[1,2,3,4,5,5+i_P]]

        # Solve and summarise PDE
        sol = SolvePDE(D,λ,K,α,Tfine,L=L,t₀=MLE_Combined_Settings[:t₀],u₀=u₀)
        Xplot = zeros(size(Tfine))
        for (i,t) ∈ enumerate(Tfine)
            Xplot[i] = SummaryStatistics(sol(t),L=L,τₐ=τ*K)[Sᵢ]
        end
        ax.plot(Tfine,Xplot,"-",color="#333333",lw=1)


    # Plot solution for each set of summary statistics
    for (i_C,Cᵢ) ∈ enumerate(MLE_Individual_Settings[:C])

        # MLE
        D,λ,K,α,τ,u₀ = GetΘ(MLE_Individual[i_C][i_P][:θ̂],MLE_Individual_Settings[:β],MLE_Individual_Settings[:i_β])

        # Solve and summarise PDE
        sol = SolvePDE(D,λ,K,α,Tfine,L=L,t₀=MLE_Individual_Settings[:t₀],u₀=u₀)
        Xplot = zeros(size(Tfine))
        for (i,t) ∈ enumerate(Tfine)
            Xplot[i] = SummaryStatistics(sol(t),L=L,τₐ=τ*K)[Sᵢ]
        end
        ax.plot(Tfine,Xplot,LineStyles[i_C],color=PoreSizeCols[i_P],lw=1)

    end

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

NumberPlots!(axs)
plt.tight_layout(pad=0.0,h_pad=1.0,w_pad=1.0)
savefig("Figures/Saved/Fig-4-Fits.pdf")
display(fig4)
