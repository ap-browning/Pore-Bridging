#=

    QuantitativeResults.jl

=#

using PoresIdentifiability
using DataFrames, JLD2, StatsBase

# Number bridged for each pore size
FinalCoverage = [Data[(Data.Day .== 28) .& (Data.PoreSize .== L),:S2_Coverage] for L ∈ P]
ProportionBridged = [sum(FinalCoverage[i] .> 0.9999) / length(FinalCoverage[i]) for i = 1:4]

# Mean proportion of edge density
Data[:,:ProportionOfEdge] = Data[:,:S1_Density] ./ Data[:,:S4_EdgeDensity]
ProportionEdge = [mean(Data[(Data.Day .== 28) .& (Data.PoreSize .== L),:ProportionOfEdge]) for L ∈ P]

# Average summary statistics from day 7 to 10
1 .- [mean(Data[(Data.Day .== 10) .& (Data.PoreSize .== L),:S1_Density]) / mean(Data[(Data.Day .== 7) .& (Data.PoreSize .== L),:S1_Density]) for L ∈ P]

# MLE table
@load "Results/MLE_Individual.jld2"
@load "Results/MLE_Combined.jld2"
@load "Results/Profiles_Individual.jld2"
@load "Results/Profiles_Combined.jld2"

ConfInts = DataFrame(L = Float64[], ψ = Symbol[], ψ̂ = Float64[], lwr = Float64[], upr = Float64[])
for (i_ψ,ψ) ∈ enumerate(Symbol.(["D","λ","K"]))

    # Per pore
    for (i_P,L) ∈ enumerate(P)

        # MLE
        ψ̂ = MLE_Individual[1][i_P][:θ̂][i_ψ]

        # Confidence interval
        lwr,upr = ConfidenceInterval(Profiles_Individual[1][i_P][:Ψ][i_ψ],Profiles_Individual[1][i_P][:p̂][i_ψ])

        push!(ConfInts,[L,ψ,ψ̂,lwr,upr])

    end

    # Combined
    ψ̂ = MLE_Combined[:θ̂][i_ψ]

    # Combined CI
    lwr,upr = ConfidenceInterval(Profiles_Combined[:Ψ][i_ψ],Profiles_Combined[:p̂][i_ψ])

    push!(ConfInts,[0.0,ψ,ψ̂,lwr,upr])

end

ConfIntsRounded = copy(ConfInts)
ConfIntsRounded[:,[:ψ̂,:lwr,:upr]] = round.(ConfIntsRounded[:,[:ψ̂,:lwr,:upr]],sigdigits=3)

## (Temp) export to latex table
file = open("ConfInts.tex","w")
for (i_P,L) ∈ enumerate([P;0])
    IntL = Int(L)
    if L != 0
        println(file,"\\SI{$IntL}{\\micro\\metre}")
    else
        println(file,"All")
    end
    for (i_ψ,ψ) ∈ enumerate(Symbol.(["D","λ","K","u₀"]))
        if i_ψ == 4
            continue
        end
        row = ConfInts[(ConfInts.L .== L) .& (ConfInts.ψ .== ψ),[:ψ̂,:lwr,:upr]]
        ψ̂   = round.(row.ψ̂[1],sigdigits=3)
        lwr = round.(row.lwr[1],sigdigits=3)
        upr = round.(row.upr[1],sigdigits=3)
        if i_ψ == 1
            ψ̂ = ψ̂ >= 100 ? Int(ψ̂) : ψ̂
            lwr = lwr >= 100 ? Int(lwr) : lwr
            upr = upr >= 100 ? Int(upr) : upr
        end
        display(ψ̂)
        println(file,"\t & $ψ̂ \t & ($lwr,$upr)")
    end
    println(file,"\\\\")
end
close(file)
