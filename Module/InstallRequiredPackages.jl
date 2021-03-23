#=
#
#   InstallRequiredPackages.jl
#
#   Run to install all required julia packages called by Pores
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

using Pkg;

Pkg.add("PyPlot")
Pkg.add("DifferentialEquations")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Statistics")
Pkg.add("StatsBase")
Pkg.add("LinearAlgebra")
Pkg.add("JLD2")
Pkg.add("NLopt")
Pkg.add("TimerOutputs")
