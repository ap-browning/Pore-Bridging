#=
#
#   SummaryStatistics.jl
#
#   Compute summary statistics given discretised quarter-domain PDE solution U (an N × N matrix)
#
#   Full details are available in the supplementary material document.
#
#   Alexander P. Browning
#       School of Mathematical Sciences
#       Queensland University of Technology
#       ap.browning@qut.edu.au  (institution)
#       ap.browning@icloud.com  (persistent)
#       https://alexbrowning.me
#
=#

function SummaryStatistics(U; L::Float64, τₐ::Float64=0.001)

    # Mesh size
    N, = size(U)
    Δx = L / (2 * (N - 1))

    # Loop through all "squares" (cornered by 4 nodes)
    numb = 0.0          # Cumulative cell count
    area = 0.0          # Cumulative area
    peri = 0.0          # Cumulative perimeter
    Ut   = U  .- τₐ     # +ve = visible, -ve = invisible
    Up   = Ut .> 0.0    # Visible cell material

    for i = 1:N-1
        for j = 1:N-1

            # Density on corners of the square
            V   = Ut[i:i+1,j:j+1]
            Vi  = Up[i:i+1,j:j+1]

            # Area and perimeter contribution
            acont,pcont = AreaPerimContribution(V,Vi,i == 1,j == 1)
            area += acont
            peri += pcont

            # Integrate number using quadrature (average in each square = trapezoid rule)
            numb += sum(U[i:i+1,j:j+1]) / 4

        end
    end

    # Multiple by four (bring up to whole domain) and scale by Δx
    area *= 4 * Δx^2
    peri *= 4 * Δx
    numb *= 4 * Δx^2

    # Return dictionary of summary statistics
    return Dict(
        :S1_Density         => numb / L^2,
        :S2_Coverage        => area / L^2,
        :S3_Circularity     => area / L^2 > 0.99 ? 1.0 : max(0.0,4π*(L^2 - area)/(peri^2 * (1 - 4π/16)) + 1 - 1 / (1 - 4π/16)),
        :S4_EdgeDensity     => U[1,1]
        )

end

# Calculate area and perimeter of unit square
function AreaPerimContribution(V,Vi,iedge,jedge)

    # Whioch case? Number of positive corners
    Npos = sum(Vi)

    # Npos == 0
    if Npos == 0

        if iedge && jedge
            acont,pcont = 0.0,2.0
        elseif iedge || jedge
            acont,pcont = 0.0,1.0
        else
            acont,pcont = 0.0,0.0
        end

    # Npos == 4
    elseif Npos == 4

        acont,pcont = 1.0,0.0

    # Npos == 3
    elseif Npos == 3

        acont,pcont = AreaPerimContribution(-V,(!).(Vi),false,false)
        acont = 1 - acont

    # Npos == 1 or 2
    else

        # Rotate to the "standard" view
        a,b,c,d = RotatePermutation(V,Vi,Npos)

        # Npos == 1
        if Npos == 1

            # Geometry
            b̃ = a / (a - b)
            c̃ = a / (a - c)

            # Compute perimeter contribution
            pcont = sqrt(b̃^2 + c̃^2)

            # Compute area contribution
            acont = 0.5 * b̃ * c̃

            # Corner
            if iedge && jedge
                pcont += (1 - b̃) + (1 - c̃)
            elseif iedge
                pcont += (1 - b̃)
            elseif jedge
                pcont += (1 - c̃)
            end

        # Npos == 2
        else

            # Geometry
            ã = a / (a - b)
            c̃ = c / (c - d)

            # Compute perimeter contribution
            pcont = sqrt(1 + (c̃ - ã)^2)

            # Compute area contribution
            acont = ã + (c̃ - ã) / 2

        end

    end

    return acont,pcont

end

# Rotation permutations (i.e., identical configerations, just rotated)
function RotatePermutation(V,Vi,Npos)

    if Npos == 2
        if Vi == [1 0; 1 0]
            a,b,c,d = V[[2,4,1,3]]
        elseif Vi == [0 0; 1 1]
            a,b,c,d = V[[4,3,2,1]]
        elseif Vi == [0 1; 0 1]
            a,b,c,d = V[[3,1,4,2]]
        else
            a,b,c,d = V
        end
    else
        if Vi == [0 0; 1 0]
            a,b,c,d = V[[2,4,1,3]]
        elseif Vi == [0 0; 0 1]
            a,b,c,d = V[[4,3,2,1]]
        elseif Vi == [0 1; 0 0]
            a,b,c,d = V[[3,1,4,2]]
        else
            a,b,c,d = V
        end
    end

    return a,b,c,d

end
