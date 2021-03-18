
using PoresIdentifiability
using PyPlot, StatsBase, JLD2
include("../Preferences.jl")

# Create figure
figS3,axs = subplots(1,3,figsize=(6.7,2.2))

# Load results
@load "Results/StandardDeviations.jld2"

# Axis limits to use
ylimits = [[0.0,0.0015],[0.0,0.5],[0.0,0.6]]
ylimits = [diff(x)[1]*[-0.05,0.05] + x for x ∈ ylimits]

for (i,Sᵢ) in enumerate(S[1:3])

    # Plot data
    axs[i].plot(μσ_df[Sᵢ].μ,μσ_df[Sᵢ].σ,"k.")

    # Plot fit
    xplot = range(Stks[Sᵢ][1],Stks[Sᵢ][end],length=200)
    axs[i].plot(xplot,σ[Sᵢ].(xplot),"r:")

    axs[i].set_xlim(Slim[Sᵢ])
    axs[i].set_ylim(ylimits[i])
    axs[i].set_xticks(Stks[Sᵢ])
    if i == 1
        axs[i].set_ylabel("σ")
    end
    axs[i].set_xlabel("μ")
    axs[i].set_title(Snms[Sᵢ])

end

NumberPlots!(axs)
plt.tight_layout()
savefig("Figures/Supplementary/Fig-S3-StandardDeviations.pdf")
display(figS3)
