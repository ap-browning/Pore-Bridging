#=
#
#   ProduceMainResults.jl
#
#   Produce all results for main paper.
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

# Note 1: The runtime is significant (multiple days, see README for more information)
# Note 2: Results should be run in this order (the profiles script, for example, depends on the MLE results)

# 1. Standard deviations
@time include("StandardDeviations.jl")

# 2. Maximum likelihood estimates
@time include("MLE.jl")

# 3. Profiles (univariate)
@time include("Profiles.jl")

# 4. Profiles (bivariate)
@time include("Profiles2D.jl")
