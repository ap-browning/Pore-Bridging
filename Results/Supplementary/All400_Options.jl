#=
#
#   All400_Options.jl
#
#   Settings for MLE analysis using all 400 μm data
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

## Fixed parameters
#      α    τ
β   = [1.0, 0.5 ]
i_β = [4,   5   ]

## Parameter bounds
#      D        λ       K       u₀
lb  = [10.0,    0.01,   0.002,  1e-5    ]
ub  = [2000.0,  2.0,    0.005,  0.002   ]

## Initial guess for MLE
#  Note: Since we first use a global optimisation algorithm,
#        this choice shouldn't matter
θ₀  = lb + rand(length(lb)) .* (ub - lb)

## Start time for PDE
t₀  = 4.0

## Tolerance to use in local maximisation algorithms
ftol_abs = 1e-4

## Runtime for global maximisation algorithms
maxtime  = 6.0 * 3600.0     # 6 hours

## Include all 400 μm data here
Data₁   = copy(Data)

## Profiles resolution
profile_res_uni = 200   # For univariate profiles
profile_res_biv = 50    # For bivariate profiles
