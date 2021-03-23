## PRODUCE SUPPLEMENTARY RESULTS
#
# @time include("All400_MLE.jl")
#
@time include("Fisher_MLE.jl")
#
@time include("Circularity_MLE.jl")

#@time include("All400_Profiles.jl")

@time include("Fisher_Profiles.jl")
