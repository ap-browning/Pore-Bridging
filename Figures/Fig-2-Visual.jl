
using PoresIdentifiability
using JLD2, PyPlot, StatsBase
include("Preferences.jl")

# Times to plot
Tplot = collect(4.0:3.0:28.0)

# Create figure
fig2,axs = subplots(4,length(Tplot),figsize=(6.7,3))
figC,axC = subplots(1,1)    # For colorbar

# Load results
@load "Results/MLE_Individual.jld2"

# Store normalised density in centre (fix plotting artifacts)
Ucentre = zeros(4,length(Tplot))

# Loop over P
for (i_P,L) ∈ enumerate(P)


    # MLE
    D,λ,K,α,τ,u₀ = GetΘ(MLE_Individual[1][i_P][:θ̂],MLE_Individual_Settings[:β],MLE_Individual_Settings[:i_β])

    # Solve PDE
    sol = SolvePDE(D,λ,K,α,Tplot,L=L,t₀=MLE_Individual_Settings[:t₀],u₀=u₀)

    # Loop and plot through time
    for (i_T,t) ∈ enumerate(Tplot)

        # Current axis
        ax    = axs[i_P,i_T]

        # Get full solution
        U_q   = sol(t)
        U     = [U_q U_q[:,end:-1:1]; U_q[end:-1:1,:] U_q[end:-1:1,end:-1:1]]
        Û     = min.(U ./ K, 1.0)

        # Plot
        c = ax.contourf(Û,levels=0.0:0.1:1.0,vmin=0.0,vmax=2.0,cmap="binary",extend="min")

        # Plot tissue boundary (if required)
        if (maximum(Û) > τ) & (U_q[end,end] / K < τ)
            ax.contour(Û,levels=[τ],colors=[PoreSizeCols[i_P]],linewidths=[1.5])
        end

        # Style axis
        ax.set_aspect("equal","box")
        ax.set_xticks([])
        ax.set_yticks([])
        if i_T == 1
            IntL = Int(L)
            ax.set_ylabel("$IntL μm")
        end
        if i_P == 1
            Intt = Int(t)
            ax.set_title("$Intt d")
        end

        # Plot colorbar
        if i_P == 1 & i_T == 1
            colorbar(c,axC)
        end

    end


end

plt.tight_layout(pad=1.0,h_pad=0.5,w_pad=-2.0)
figure(fig2.number)
savefig("Figures/Fig-2-Visual.pdf")
display(fig2)

# Plot to get colorbar
figure(figC.number)
savefig("Figures/Fig-2-Visual-Colorbar.pdf")
