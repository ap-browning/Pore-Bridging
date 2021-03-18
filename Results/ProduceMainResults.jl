## PRODUCE MAIN RESULTS

# Standard deviations
@time include("StandardDeviations.jl")

# Maximum likelihood estimates
@time include("MLE.jl")

# Profiles (bivariate)
@time include("Profiles2D.jl")

# Profiles (univariate)
@time include("Profiles.jl")
