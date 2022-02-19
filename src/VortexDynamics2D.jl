module VortexDynamics2D

using LinearAlgebra



"""
    biotsavart(px, py, Γ; aggregate=sum)

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
function biotsavart(px::Vector{Float64}, py::Vector{Float64}, Γ::Vector{Float64}; aggregate=sum, core=1.0e-6)

    ppx = @view px[:,:]
    ppy = @view py[:,:]
    
    # represent position vectors as row/column tensors
    ppx = reshape(ppx, :, 2)
    ppy = reshape(ppy, :, 2)

    rx = permutedims( ppx[:,:,:], [1,3,2])
    ry = permutedims( ppy[:,:,:], [3,1,2])


    D = rx .- ry # tensorized difference between all dimensions
    L = sum( abs2, D, dims=3) .+ 0.0 # compute euclidean norm for each D(i,j,:)

    # regularize the vortices by capping the stream function inside the core
    # (this would result in zero velocity there --- this ensures that the vortex
    # isn't moving itself)
    L[L .< core] .= core

    # compute the Biot--Savart
    #Ψ = -log.(L).*Γ/2/π

    gradΨ  =  cat( -D[:,:,2], D[:,:,1]; dims=3) ./ (L.^2) .* Γ /2/π; 

    V = sum( gradΨ, dims=2)
    V = permutedims( V[:,:,:], [1,3,2] )
    V = reshape(V, :, 1)


    

end


end
