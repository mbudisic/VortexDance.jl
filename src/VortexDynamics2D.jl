module VortexDynamics2D

using LinearAlgebra
using Infiltrator


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
function biotsavart(pv, pe, Γ; aggregate=:sum, core=1.0e-6)

    pvx = @view pv[1:2:end]
    pvy = @view pv[2:2:end]

    pex = @view pe[1:2:end]
    pey = @view pe[2:2:end]

    # ensure that number of vertices is the same as number of circulations
    @assert size(pvx) == size(Γ)

    D = -cat( pvy .- transpose( pey ), -( pvx .- transpose( pex ) ); dims=3  )

    L = dropdims( sum( abs2, D, dims=3) .+ 0.0; dims=3) # compute euclidean norm for each D(i,j,:)
    
    strongest_vortex = findmax(L; dims = 1)

    D[ L .< core .^2, : ] .= 0.0
    L[ L .< core .^2 ] .= core .^2 
    # compute the Biot--Savart
    #Ψ = -log.(L).*Γ/2/π

    perp_gradΨ  =  - (Γ/2/π) .* D ./ L ; 
    # perp_gradΨ[sel,:] .= 0.0

    # size(perp_gradΨ) = # vortices x # evaluation points x 2 (x,y)

    # add all vortex influences  size(V) = 1 x #evaluation points x 2
    if aggregate == :sum
        V = dropdims( sum( perp_gradΨ, dims=1 ); dims=1 )
        V = transpose(V)
    elseif aggregate == :max
        V = dropdims( perp_gradΨ[strongest_vortex[2],:]; dims=1)
        V = transpose(V)
    else
        error("Unknown aggregation")
    end

    V = V[:]
    

end


end
