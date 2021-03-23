#=
#
#   Timings.jl
#
#   Create data frame of result timings for display in README file.
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

using JLD2, DataFrames

@load "Results/Saved/MLE_Individual.jld2"
@load "Results/Saved/MLE_Combined.jld2"
@load "Results/Saved/Profiles_Individual.jld2"
@load "Results/Saved/Profiles_Combined.jld2"
@load "Results/Saved/Profiles2D.jld2"

Timings = DataFrame(Script = String[], Description = String[], Hours = Float64[])

# Computing MLE
push!(Timings, (
    "MLE.jl",
    "Compute individual and combined MLEs",
    MLE_Individual_Settings[:walltime] / 3600 + 
    MLE_Combined_Settings[:walltime] / 3600
    ))

# Profiles
push!(Timings, (
    "Profiles.jl",
    "Compute individual and combined profile likelihoods",
    Profiles_Individual_Settings[:walltime] / 3600 + 
    Profiles_Combined_Settings[:walltime] / 3600
    ))

# Bivariate profiles
push!(Timings, (
    "Profiles2D.jl",
    "Compute bivariate profiles",
    Profiles2D_Settings[:walltime] / 3600
    ))