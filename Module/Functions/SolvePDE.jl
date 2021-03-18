#=
    PorePDE.jl
    Solve Porous-Fisher equation on an L × L domain
    author:  Alexander P Browning
    contact: ap.browning@icloud.com
=#

function SolvePDE(D::Float64,λ::Float64,K::Float64,α::Float64,T::Array{Float64,1};
                 L::Float64=300.0,t₀::Float64=0.0,u₀::Float64=0.001,N::Int=51,method=Tsit5())

    # Calculate constants

        # Mesh size
        Δx = L / (2N - 2)

        # β
        β  = D / (2 * K^α * Δx^2)

        # Boundary lookup functions (reflecting)
        IndexUp = i -> i != N ? i + 1 : N - 1
        IndexDn = i -> i != 1 ? i - 1 : 2

    # Discretised PDE
    function f!(dudt::Array{Float64,2},U::Array{Float64,2},p,t)

        if α == 0.0
            # Fisher
            V = ones(size(U))
        elseif α == 1.0
            # Porous-Fisher
            V = U
        else
            # Generalised Fisher
            V = @. abs(U)^α
        end


        # Logistic growth
        @. dudt = λ * U * (1.0 - U / K)

        # Inner nodes
        @inbounds @simd for i = 2:N-1
            for j = 2:N-1
                dudt[i,j] +=
                            β * ((V[i,j] + V[i+1,j]) * (U[i+1,j] - U[i,j]) - (V[i-1,j] + V[i,j]) * (U[i,j] - U[i-1,j])) +
                            β * ((V[i,j] + V[i,j+1]) * (U[i,j+1] - U[i,j]) - (V[i,j-1] + V[i,j]) * (U[i,j] - U[i,j-1]))
            end
        end

        # Boundary nodes
        @inbounds @simd for i = [N]
            for j = 2:N-1
                dudt[i,j] +=
                            β * ((V[i,j] + V[IndexUp(i),j]) * (U[IndexUp(i),j] - U[i,j]) - (V[IndexDn(i),j] + V[i,j]) * (U[i,j] - U[IndexDn(i),j])) +
                            β * ((V[i,j] + V[i,IndexUp(j)]) * (U[i,IndexUp(j)] - U[i,j]) - (V[i,IndexDn(j)] + V[i,j]) * (U[i,j] - U[i,IndexDn(j)]))
            end
        end

        # Boundary nodes
        @inbounds @simd for i = 2:N
            for j = [N]
                dudt[i,j] +=
                            β * ((V[i,j] + V[IndexUp(i),j]) * (U[IndexUp(i),j] - U[i,j]) - (V[IndexDn(i),j] + V[i,j]) * (U[i,j] - U[IndexDn(i),j])) +
                            β * ((V[i,j] + V[i,IndexUp(j)]) * (U[i,IndexUp(j)] - U[i,j]) - (V[i,IndexDn(j)] + V[i,j]) * (U[i,j] - U[i,IndexDn(j)]))
            end
        end

    end

    # Fill initial condition
    U::Array{Float64,2} = zeros(N,N)
    U[1,:] .= u₀
    U[:,1] .= u₀

    # Solve
    return solve(ODEProblem(f!,U,(t₀,T[end])),method,saveat=T)

end
