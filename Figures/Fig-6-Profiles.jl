
using PoresIdentifiability
using PyPlot, StatsBase, JLD2
include("Preferences.jl")

# Create figure
fig6,axs = subplots(5,4,figsize=(6.7,6.0))

# Load results
@load "Results/Profiles_Individual.jld2"
@load "Results/Profiles_Combined.jld2"
@load "Results/MLE_Individual.jld2"
@load "Results/MLE_Combined.jld2"

# Options
LineStyles = ["-","--"]
Ψnames     = ["D","λ","K","u₀"]
Ψticks     = [0.:500.:2000., 0.:0.5:2., 0.002:0.001:0.005, 0.0:0.0005:0.002]

# Loop over P and parameters
for (i_P,L) ∈ enumerate([P;0]), i_ψ = 1:4

    # Current axis
    ax = axs[i_P,i_ψ]

    # Plot profile for each set of summary statistics
    for (i_C,Cᵢ) ∈ enumerate(Profiles_Individual_Settings[:C])

        if L == 0.0
            continue
        end

        # Profile
        Ψᵢ = Profiles_Individual[i_C][i_P][:Ψ][i_ψ]
        p̂ᵢ = Profiles_Individual[i_C][i_P][:p̂][i_ψ]
        p̂ᵢ .-= maximum(p̂ᵢ)
        p̂ᵢ = max.(p̂ᵢ,-100)

        # Plot
        ax.fill_between(Ψᵢ,-6,p̂ᵢ,color=PoreSizeCols[i_P],alpha=(i_C==2 ? 0.2 : 0.1))
        ax.plot(Ψᵢ,p̂ᵢ,LineStyles[i_C],c=PoreSizeCols[i_P])

        if i_C == 1

            θ̂ = MLE_Individual[i_C][i_P][:θ̂][i_ψ]
            ax.plot(θ̂,0.0,".",color=PoreSizeCols[i_P])

        end

    end

    # Plot combined (coded as L = 0)
    if L == 0 && i_ψ <= 3

        # Profile
        Ψᵢ = Profiles_Combined[:Ψ][i_ψ]
        p̂ᵢ = Profiles_Combined[:p̂][i_ψ]
        p̂ᵢ .-= maximum(p̂ᵢ)
        p̂ᵢ = max.(p̂ᵢ,-100)

        # Plot
        ax.fill_between(Ψᵢ,-6,p̂ᵢ,color="#555555",alpha=0.1)
        ax.plot(Ψᵢ,p̂ᵢ,LineStyles[1],c="#555555")

        θ̂ = MLE_Combined[:θ̂][i_ψ]
        ax.plot(θ̂,0.0,".",color="#555555")

    end

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
    if i_P == 5
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
savefig("Figures/Fig-6-Profiles.pdf")
display(fig6)
