
using Pores
using PyPlot, StatsBase, JLD2
include("../Preferences.jl")

# Create figure
figS4a,axs = subplots(1,4,figsize=(6.7,1.8))

# Load results
@load "Results/Supplementary/Saved/All400_Profiles_Individual.jld2"
@load "Results/Supplementary/Saved/All400_MLE_Individual.jld2"

# Options
LineStyles = ["-","--"]
Ψnames     = ["D","λ","K","u₀"]
Ψticks     = [0.:500.:2000., 0.:0.5:2., 0.002:0.001:0.005, 0.0:0.0005:0.002]

# L = 400
i_P,L = 2,400.0
i_C   = 1

# Loop over P and parameters
for i_ψ = 1:4

    # Current axis
    ax = axs[i_ψ]

        # Profile
        Ψᵢ = All400_Profiles_Individual[i_C][i_P][:Ψ][i_ψ]
        p̂ᵢ = All400_Profiles_Individual[i_C][i_P][:p̂][i_ψ]
        p̂ᵢ .-= All400_MLE_Individual[i_C][i_P][:L̂]
        p̂ᵢ = max.(p̂ᵢ,-100)

        # Plot
        ax.fill_between(Ψᵢ,-6,p̂ᵢ,color=PoreSizeCols[i_P],alpha=0.3)
        ax.plot(Ψᵢ,p̂ᵢ,LineStyles[i_C],c=PoreSizeCols[i_P])

        θ̂ = All400_MLE_Individual[i_C][i_P][:θ̂][i_ψ]
        ax.plot(θ̂,0.0,".",color=PoreSizeCols[i_P])

    # 95% confidence interval cutoff
    xlimit = ax.get_xlim()
    ax.plot(xlimit,[-1.92,-1.92],"k:")
    ax.set_xlim(xlimit)

    # Style axes
    ax.set_ylim([-4.2,0.2])
    ax.set_xticks(Ψticks[i_ψ])
    if i_ψ == 1
        IntL = Int(L)
        ax.set_ylabel("$IntL μm")
    else
        ax.set_yticklabels([])
    end
    ax.set_xlabel(Ψnames[i_ψ])

    if i_ψ == 4
        ax.set_xlim([-0.0001,0.0021])
    end

end

NumberPlots!(axs)
plt.tight_layout(pad=1.0)
savefig("Figures/Supplementary/Saved/Fig-S4b-All400_Profiles.pdf")
display(figS4a)
