#=

    PoresIdentifiability.jl

    A Julia module for the pore identifiability project

    author:  Alexander P Browning
    contact: ap.browning@icloud.com
    push!(LOAD_PATH,"C:/Users/browniap/OneDrive - Queensland University of Technology/Desktop/Pores/Module")

=#
module PoresIdentifiability

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

    include("Functions/SolvePDE.jl")
    include("Functions/SummaryStatistics.jl")
    include("Functions/Maximise.jl")
    include("Functions/Profile.jl")
    include("Functions/LogLikelihood.jl")

    export SolvePDE, SummaryStatistics, Maximise, FitModel, Profile, LogLikelihood, GetΘ, ConfidenceInterval

    # Import data
    Data = DataFrame(CSV.File("Module/Data.csv"))

    # Set of all pore sizes:            L ∈ P
    P = Float64.(sort(unique(Data.PoreSize)))

    # Set of all observation times:     t ∈ T
    T = Float64.(sort(unique(Data.Day)))

    # Set of all summary statistics:    Sᵢ ∈ S
    S = [:S1_Density,:S2_Coverage,:S3_Circularity,:S4_EdgeDensity]

    export Data, P, T, S

end
