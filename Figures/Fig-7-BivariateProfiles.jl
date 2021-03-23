
using Pores
using JLD2, PyPlot, StatsBase
include("Preferences.jl")

# Create figure
fig7,axs = subplots(2,4,figsize=(6.7,3.0))
fcbar,acbar = subplots(1,1) # For colorbar

# Load results
@load "Results/Saved/MLE_Individual.jld2"
@load "Results/Saved/Profiles2D.jld2"

# Loop over P
for (i_P,L) ∈ enumerate(P)

     Ψ = Profiles2D[i_P][:Ψ]
     θ̂ = MLE_Individual[1][i_P][:θ̂]

     # Plot 1: λ vs D
     bv_p̂₁ = copy(Profiles2D[i_P][:bv_p̂₁])
     bv_p̂₁ .-= maximum(bv_p̂₁)
     bv_p̂₁ = max.(-11.0,bv_p̂₁)

     ax = axs[1,i_P]
     ax.contourf(Ψ[2],Ψ[1],bv_p̂₁,levels=-11:0.5:0)
     ax.contour(Ψ[2],Ψ[1],bv_p̂₁,levels=[-3.00],colors=["#ffffff"])

     ax.plot(θ̂[2],θ̂[1],"rx")

     # Style axes
     if i_P == 1
         ax.set_ylabel("D")
     else
         ax.set_yticklabels([])
     end
     ax.set_title("$L μm")
     ax.set_xlim([0.0,2.0])
     ax.set_ylim([0.0,2000.0])
     ax.set_xticklabels([])
     ax.set_aspect(diff([ax.get_xlim()...])/diff([ax.get_ylim()...]))

     # Plot 2: λ vs K
     bv_p̂₂ = copy(Profiles2D[i_P][:bv_p̂₂])
     bv_p̂₂ .-= maximum(bv_p̂₂)
     bv_p̂₂ = max.(-11.0,bv_p̂₂)

     ax = axs[2,i_P]
     c2 = ax.contourf(Ψ[2],Ψ[3],bv_p̂₂',levels=-11:1:0)
     ax.contour(Ψ[2],Ψ[3],bv_p̂₂',levels=[-6.00],colors=["w"])

     ax.plot(θ̂[2],θ̂[3],"rx")

     # Style axes
     if i_P == 1
         ax.set_ylabel("K")
     else
         ax.set_yticklabels([])
     end
     ax.set_xlim([0.0,2.0])
     ax.set_ylim([0.002,0.005])
     ax.set_xlabel("λ")
     ax.set_aspect(diff([ax.get_xlim()...])/diff([ax.get_ylim()...]))


     # Plot colorbar
     if i_P == 1
         colorbar(c2,acbar)
     end

end

NumberPlots!(axs)
figure(fig7.number);
plt.tight_layout(pad=0.0,h_pad=0.0,w_pad=-10.0)
savefig("Figures/Saved/Fig-7-BivariateProfiles.pdf")
figure(fcbar.number); savefig("Figures/Saved/Fig-7-BivariateProfiles-Colorbar.pdf")
display(fig7)
