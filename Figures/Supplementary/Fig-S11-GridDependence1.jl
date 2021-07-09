
using Pores
using JLD2, PyPlot, StatsBase
include("../Preferences.jl")

# Times to plot (per initial condition)
Tplot = [14.0,28.0]

# Create figure
figS11,axs = subplots(3,4,figsize=(6.7,5))

# Load results
@load "Results/Saved/MLE_Combined.jld2"

# Second grid spacing to try
N₂ = 151

# Loop over P
for (i_P,L) ∈ enumerate(P)

    # MLE
    D,λ,K,α,τ,u₀ = GetΘ(MLE_Combined[:θ̂],MLE_Combined_Settings[:β],MLE_Combined_Settings[:i_β])[[1,2,3,4,5,5+i_P]]

    # Solve PDE (N = 51)
    sol1 = SolvePDE(D,λ,K,α,Tplot,L=L,t₀=MLE_Combined_Settings[:t₀],u₀=u₀)
    
    # Solve PDE (N = N₂)
    sol2 = SolvePDE(D,λ,K,α,Tplot,L=L,t₀=MLE_Combined_Settings[:t₀],u₀=u₀,N=N₂)

    # Plot density profiles
    for t ∈ Tplot

        U1   = sol1(t)[end,:]
        X1   = range(0.0,L / 2,length=length(U1))

        U2   = sol2(t)[end,:]
        X2   = range(0.0,L / 2,length=length(U2))

        # Plot
        axs[1,i_P].plot(X1,U1,color=PoreSizeCols[i_P])
        axs[1,i_P].plot(X2,U2,":k")

    end

    axs[1,i_P].set_ylim([-0.00025,0.00425])
    axs[1,i_P].set_ylabel("Density")
    axs[1,i_P].set_title("$(Int(L)) μm")
    axs[1,i_P].set_xticks(range(0.0,L/2,length=3))
    axs[1,i_P].set_xlabel("x (μm)")

    # Plot solution for N = 51 at t = t₁
    U   = sol1(Tplot[1])[end:-1:1,:]
    Û     = min.(U ./ K, 1.0)
    axs[2,i_P].contourf(Û,levels=0.0:0.1:1.0,vmin=0.0,vmax=2.0,cmap="binary",extend="min")
    if (maximum(Û) > τ) & (U[1,end] / K < τ)
        axs[2,i_P].contour(Û,levels=[τ],colors=[PoreSizeCols[i_P]],linewidths=[1.5])
    end
    axs[2,i_P].set_aspect("equal","box")
    axs[2,i_P].set_xticks([])
    axs[2,i_P].set_yticks([])

    # Plot solution for N = N₂ at t = t₁
    U   = sol2(Tplot[1])[end:-1:1,:]
    Û     = min.(U ./ K, 1.0)
    axs[3,i_P].contourf(Û,levels=0.0:0.1:1.0,vmin=0.0,vmax=2.0,cmap="binary",extend="min")
    if (maximum(Û) > τ) & (U[1,end] / K < τ)
        axs[3,i_P].contour(Û,levels=[τ],colors=[PoreSizeCols[i_P]],linewidths=[1.5])
    end
    axs[3,i_P].set_aspect("equal","box")
    axs[3,i_P].set_xticks([])
    axs[3,i_P].set_yticks([])

end

NumberPlots!(axs)
plt.tight_layout(pad=0.0,h_pad=1.0,w_pad=1.0)
display(figS11)

savefig("Figures/Supplementary/Fig-S11-GridDependence1.pdf")