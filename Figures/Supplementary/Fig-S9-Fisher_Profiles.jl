
using Pores
using PyPlot, StatsBase, JLD2
include("../Preferences.jl")

# Create figure
fig6,axs = subplots(4,4,figsize=(6.7,5.0))

# Load results
@load "Results/Supplementary/Saved/Fisher_Profiles_Individual.jld2"
@load "Results/Supplementary/Saved/Fisher_MLE_Individual.jld2"

# Options
LineStyles = ["-","--"]
Ψnames     = ["D","λ","K","u₀"]
Ψticks     = [0.:100.:500., 0.:0.5:2., 0.002:0.001:0.005, 0.0:0.0005:0.002]

# Loop over P and parameters
for (i_P,L) ∈ enumerate(P), i_ψ = 1:4

    # Current axis
    ax = axs[i_P,i_ψ]

    # Plot profile for each set of summary statistics
    for (i_C,Cᵢ) ∈ enumerate(Fisher_Profiles_Individual_Settings[:C][[1]])

        # Profile
        Ψᵢ = Fisher_Profiles_Individual[i_C][i_P][:Ψ][i_ψ]
        p̂ᵢ = Fisher_Profiles_Individual[i_C][i_P][:p̂][i_ψ]
        p̂ᵢ .-= Fisher_MLE_Individual[i_C][i_P][:L̂]
        p̂ᵢ = max.(p̂ᵢ,-100)

        # Plot
        ax.fill_between(Ψᵢ,-6,p̂ᵢ,color=PoreSizeCols[i_P],alpha=0.3)
        ax.plot(Ψᵢ,p̂ᵢ,LineStyles[i_C],c=PoreSizeCols[i_P])

        if i_C == 1

            θ̂ = Fisher_MLE_Individual[i_C][i_P][:θ̂][i_ψ]
            ax.plot(θ̂,0.0,".",color=PoreSizeCols[i_P])

        end

    end

    # 95% confidence interval cutoff
    xlimit = ax.get_xlim()
    ax.plot(xlimit,[-1.92,-1.92],"k:")
    ax.set_xlim(xlimit)

    # Adjust D limit for Fisher (crop to (0,500))
    if i_ψ == 1
        ax.set_xlim([-25,525])
    end

    # Style axes
    ax.set_ylim([-4.2,0.2])
    ax.set_xticks(Ψticks[i_ψ])
    if i_ψ == 1
        IntL = Int(L)
        ax.set_ylabel("$IntL μm")
    else
        ax.set_yticklabels([])
    end
    if i_P == 4
        ax.set_xlabel(Ψnames[i_ψ])
    else
        ax.set_xticklabels([])
    end

    if i_ψ == 4
        ax.set_xlim([-0.0001,0.0021])
    end

end

NumberPlots!(axs)
plt.tight_layout(pad=0.0,h_pad=0.2,w_pad=0.2)
savefig("Figures/Supplementary/Saved/Fig-S9-Fisher_Profiles.pdf")
display(fig6)
