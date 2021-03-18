#=
#
#   Maximise.jl
#
#   Perform bounded maximisation using methods in NLopt.jl
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au
#       https://alexbrowning.me
#
=#
@doc """
    Maximise(fun::Function,θ0::Array{Float64,1};lb::Array{Float64,1}=zeros(length(θ0)),ub::Array{Float64,1}=fill(Inf,length(θ0)),f_tol::Float64=1e-2,method=:LN_COBYLA)

Maximises function `fun` (of the single vector variable `θ0`) within bounds
    specified by `lb` (lower) and `ub` (upper). Uses `NLopt.jl` routine given by
    `method` where the tolerence is given by `f_tol`.

Inputs:\n
    `fun`       - Function to maximise, of variable `θ`, fun(θ) : R^N → R
    `θ0`        - Initial guess
    `lb`        - lower bounds (default: zeros)
    `ub`        - upper bounds (default: Inf)
    `ftol_rel`  - tolerance (default: 0.0)
    `ftol_abs`  - tolerance (default: 0.0)
    `method`    - method to use (default: COBYLA)

Outputs:\n
    `θ̂` - argmax(fun(θ))
    `L̂` - max(fun(θ))
"""
function Maximise(fun::Function,θ0::Array{Float64,1},
                  method::Symbol=:LN_BOBYQA;
                  lb::Array{Float64,1}=zeros(length(θ0)),
                  ub::Array{Float64,1}=fill(Inf,length(θ0)),
                  maxtime::Float64=0.0,
                  ftol_rel::Float64=0.0,
                  ftol_abs::Float64=((maxtime == 0.0) & (ftol_rel == 0.0)) ? 1e-3 : 0.0,
                  kwargs...)

    while any(ub .<= θ0)
        idx = findfirst(ub .<= θ0)
        θ0[idx] = ub[idx] - 1e-9
    end
    while any(lb .>= θ0)
        idx = findfirst(lb .>= θ0)
        θ0[idx] = lb[idx] + 1e-9
    end

    # Map to bivariate function that write the derivative also (not used)
    tomax = (θ,∂θ) -> fun(θ)

    # Create Opt object, enter options
    opt = Opt(method,length(θ0))
    opt.max_objective = tomax
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.local_optimizer = Opt(:LN_BOBYQA, length(θ0))
    if ftol_abs != 0.0
        opt.ftol_abs = ftol_abs
    end
    if ftol_rel != 0.0
        opt.ftol_rel = ftol_rel
    end
    if maxtime != 0.0
        opt.maxtime = maxtime
    end

    # Maximise
    res = optimize(opt,θ0)

    # Return θ̂, L̂
    return res[[2,1]]

end

function Maximise(fun::Function,θ0::Array{Float64,1},
                  method::Array{Symbol,1};kwargs...)

    # Find global maxima
    θ̂₁, = Maximise(fun,θ0,method[1];kwargs...)

    # Find and return local maxima
    return Maximise(fun,θ̂₁,method[2];kwargs...)

end

function FitModel(fun::Function,θ0::Array{Float64,1},x,y)

    function ToMax(θ)
        nSS = 0.0
        for i = 1:length(x)
            if isfinite(y[i])
                nSS -= (y[i] - fun(x[i],θ))^2
            end
        end
        return nSS
    end

    θ, = Maximise(ToMax,θ0,lb=fill(-Inf,length(θ0)),ftol_rel=1e-6)

    return θ

end
