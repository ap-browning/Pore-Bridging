
using Pores
using JLD2, PyPlot, StatsBase
include("../Preferences.jl")

# Times to plot (per initial condition)
Tplot1 = collect(range(4.0,28.0,length=7))
Tplot2 = collect(range(4.0,28.0,length=7))
Tplot3 = collect(range(19.0,22.0,length=7))

# Create figure
figS10,axs = subplots(3,length(Tplot1),figsize=(10,4))

# Load results
@load "Results/Saved/MLE_Combined.jld2"

# Grid spacing to try
N = 151

# Pore size
i_P,L = 1,300.0

    # MLE
    D,λ,K,α,τ,u₀ = GetΘ(MLE_Combined[:θ̂],MLE_Combined_Settings[:β],MLE_Combined_Settings[:i_β])[[1,2,3,4,5,5+i_P]]

    # Solve PDE (Dirichlet)
    sol1 = SolvePDE(D,λ,K,α,Tplot1,L=L,t₀=MLE_Combined_Settings[:t₀],u₀=u₀,N=N)
    
    # Solve PDE (Neumann, set 1)
    sol2 = SolvePDE_Neumann(D,λ,K,α,Tplot2,L=L,t₀=MLE_Combined_Settings[:t₀],u₀=u₀,N=N)

    # Solve PDE (Neumann, set 2)
    sol3 = SolvePDE_Neumann(D,λ,K,α,Tplot3,L=L,t₀=MLE_Combined_Settings[:t₀],u₀=u₀,N=N)

# Loop over Tplot1
for i_T = 1:length(Tplot1)
    
    # Plot solution for Dirichlet boundary conditions
    U   = sol1(Tplot1[i_T])[end:-1:1,:]
    Û     = min.(U ./ K, 1.0)
    axs[1,i_T].contourf(Û,levels=0.0:0.1:1.0,vmin=0.0,vmax=2.0,cmap="binary",extend="min")
    if (maximum(Û) > τ) & (U[1,end] / K < τ)
        axs[1,i_T].contour(Û,levels=[τ],colors=[PoreSizeCols[i_P]],linewidths=[1.5])
    end
    axs[1,i_T].set_aspect("equal","box")
    axs[1,i_T].set_xticks([])
    axs[1,i_T].set_yticks([])
    axs[1,i_T].set_title("$(Tplot1[i_T]) d")

    # Plot solution for Neumann boundary conditions
    U   = sol2(Tplot2[i_T])[end:-1:1,:]
    Û     = min.(U ./ K, 1.0)
    axs[2,i_T].contourf(Û,levels=0.0:0.1:1.0,vmin=0.0,vmax=2.0,cmap="binary",extend="min")
    if (maximum(Û) > τ) & (U[1,end] / K < τ)
        axs[2,i_T].contour(Û,levels=[τ],colors=[PoreSizeCols[i_P]],linewidths=[1.5])
    end
    if (maximum(Û) > 0.1) & (U[1,end] / K < 0.1)
        axs[2,i_T].contour(Û,levels=[0.1],linestyles=":",colors=[PoreSizeCols[i_P]],linewidths=[1.5])
    end
    axs[2,i_T].set_aspect("equal","box")
    axs[2,i_T].set_xticks([])
    axs[2,i_T].set_yticks([])

    # Plot solution for Neumann boundary conditions
    U   = sol3(Tplot3[i_T])[end:-1:1,:]
    Û     = min.(U ./ K, 1.0)
    axs[3,i_T].contourf(Û,levels=0.0:0.1:1.0,vmin=0.0,vmax=2.0,cmap="binary",extend="min")
    if (maximum(Û) > τ) & (U[1,end] / K < τ)
        axs[3,i_T].contour(Û,levels=[τ],colors=[PoreSizeCols[i_P]],linewidths=[1.5])
    end
    if (maximum(Û) > 0.1) & (U[1,end] / K < 0.1)
        axs[3,i_T].contour(Û,levels=[0.1],linestyles=":",colors=[PoreSizeCols[i_P]],linewidths=[1.5])
    end
    axs[3,i_T].set_aspect("equal","box")
    axs[3,i_T].set_xticks([])
    axs[3,i_T].set_yticks([])
    axs[3,i_T].set_title("$(Tplot3[i_T]) d")

end

axs[1,1].set_ylabel("Dirichlet")
axs[2,1].set_ylabel("Neumann")
axs[3,1].set_ylabel("Neumann")

NumberPlots!(axs)
plt.tight_layout(pad=0.0,h_pad=1.0,w_pad=1.0)
display(figS10)

savefig("Figures/Supplementary/Saved/Fig-S10-BoundaryConditions.pdf")