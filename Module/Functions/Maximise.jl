#=
#
#   Maximise.jl
#
#   Perform bounded maximisation using methods in NLOpt.jl (https://github.com/JuliaOpt/NLopt.jl)
#
#   See https://nlopt.readthedocs.io/en/latest/ for descriptions of methods.
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

function Maximise(fun::Function,                    # Function to maximise
                  θ0::Array{Float64,1},             # Initial guess
                  method::Symbol=:LN_BOBYQA;        # Algorithm
                  lb::Array{Float64,1}=zeros(length(θ0)),       # Lower bounds
                  ub::Array{Float64,1}=fill(Inf,length(θ0)),    # Upper bounds
                  maxtime::Float64=0.0,             # Maximum walltime
                  ftol_rel::Float64=0.0,            # Relative tolerance
                  ftol_abs::Float64=((maxtime == 0.0) & (ftol_rel == 0.0)) ? 1e-3 : 0.0,    # Absolute tolerance (default)
                  kwargs...) 

    # Fix issues with some NLopt algorithms where initial guess is on the boundary
    while any(ub .<= θ0)
        idx = findfirst(ub .<= θ0)
        θ0[idx] = ub[idx] - 1e-9
    end
    while any(lb .>= θ0)
        idx = findfirst(lb .>= θ0)
        θ0[idx] = lb[idx] + 1e-9
    end

    # Map to bivariate function that write the derivative also (not used here)
    tomax = (θ,∂θ) -> fun(θ)

    # Create Opt object, enter options
    opt = Opt(method,length(θ0))

        # Objective
        opt.max_objective = tomax

        # Bounds
        opt.lower_bounds = lb       # Lower bound
        opt.upper_bounds = ub       # Upper bound

        # Some global routines require a local optimiser
        opt.local_optimizer = Opt(:LN_BOBYQA, length(θ0))

        # Figure out what tolerance to use
        if ftol_abs != 0.0
            opt.ftol_abs = ftol_abs
        end
        if ftol_rel != 0.0
            opt.ftol_rel = ftol_rel
        end
        if maxtime != 0.0
            opt.maxtime = maxtime
        end

    # Perform maximisation
    res = optimize(opt,θ0)

    # Return θ̂, L̂
    return res[[2,1]]

end

# Use two subsequent optimisation algorithms, with the same settings
function Maximise(fun::Function,θ0::Array{Float64,1},
                  method::Array{Symbol,1};kwargs...)

    # Find global maxima (using method[1])
    θ̂₁, = Maximise(fun,θ0,method[1];kwargs...)

    # Find and return local maxima (using method[1])
    return Maximise(fun,θ̂₁,method[2];kwargs...)

end

# Quickly fit model to data (i.e., fit quadratic noise function)
function FitModel(fun::Function,θ0::Array{Float64,1},x,y)

    # Maximise negative sum of squares (least-squares)
    function ToMax(θ)
        nSS = 0.0
        for i = 1:length(x)
            if isfinite(y[i])
                nSS -= (y[i] - fun(x[i],θ))^2
            end
        end
        return nSS
    end

    # Perform maximisation
    θ, = Maximise(ToMax,θ0,lb=fill(-Inf,length(θ0)),ftol_rel=1e-6)

    # Return optimal parameters
    return θ

end
