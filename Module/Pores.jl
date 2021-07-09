#=
#
#   Pores.jl
#
#   A Julia module that contains all functions and data used to produce the results.
#
#   The containing folder must be in the Julia PATH to be callable:
#       push!(LOAD_PATH,"/path/to/repo/Module")
#       using Pores
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

module Pores

    # Load packages
    using PyPlot
    using DifferentialEquations
    using CSV
    using DataFrames
    using Statistics
    using StatsBase
    using LinearAlgebra
    using .Threads
    using JLD2
    using NLopt

    # Load functions
    include("Functions/SolvePDE.jl")
    include("Functions/SolvePDE_Neumann.jl")
    include("Functions/SummaryStatistics.jl")
    include("Functions/Maximise.jl")
    include("Functions/Profile.jl")
    include("Functions/LogLikelihood.jl")

    # Functions callable once `using Pores` is called
    export SolvePDE, SolvePDE_Neumann, SummaryStatistics, Maximise, FitModel, Profile, LogLikelihood, GetΘ, ConfidenceInterval

    # Import data
    Data = DataFrame(CSV.File("Module/Data.csv"))

    # Set of all pore sizes:            L ∈ P
    P = Float64.(sort(unique(Data.PoreSize)))

    # Set of all observation times:     t ∈ T
    T = Float64.(sort(unique(Data.Day)))

    # Set of all summary statistics:    Sᵢ ∈ S
    S = [:S1_Density,:S2_Coverage,:S3_Circularity,:S4_EdgeDensity]

    # Variables callable once `using Pores` is called
    export Data, P, T, S
    
end
