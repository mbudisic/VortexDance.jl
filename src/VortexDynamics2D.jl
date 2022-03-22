module VortexDynamics2D

using LinearAlgebra



"""
    biotsavart(pv, pe, Γ; aggregate=sum)

Evaluate the stream function of a collection of 2d vortices with positions
in vector `p` and circulations in vector `Γ` using the Biot--Savart law.

`x` is the position of evaluation
The length of vector `Γ` determines the number of vortices `N`.
The length of vector `p` should be exactly `2N` with coordinates of vortices 
interlaced, i.e., `x,y,x,y,...` 

setting the keyword argument `aggregate` to something that is not `sum`
uses a different aggregation function across vortices. e.g. `aggregate=maximum` 
will use the strongest contribution, not the sum of contributions. This is nonphysical 
but can be of interest.
"""
function biotsavart(pv, pe, Γ; aggregate=sum, core=1.0e-6)

    Nv = length(pv) ÷ 2
    Ne = length(pe) ÷ 2

    L = zeros(Nv,Ne)
    D = zeros(Nv,Ne, 2)

    # perp vectors between positions
    for ix = 1:Nv
        for iy = 1:Ne
            D[ix,iy,2] = -( pv(ix)-pe(iy) )
            D[ix,iy,1] = pv(ix+Nv) - pe(iy+Ne)
        end
    end
    L = sum( abs2, D, dims=3) .+ 0.0 # compute euclidean norm for each D(i,j,:)
    L[L .< core] .= core

    # compute the Biot--Savart
    #Ψ = -log.(L).*Γ/2/π

    gradΨ  =  D ./ (L.^2) .* Γ /2/π; 

    V = sum( gradΨ, dims=2)
    V = permutedims( V[:,:,:], [1,3,2] )
    V = reshape(V, :, 1)

    V = V[:]


    

end


end
