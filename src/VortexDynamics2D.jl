module VortexDynamics2D

using LinearAlgebra
using Infiltrator

"""
    streamfunction(x, p, Γ; aggregate=sum)

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
function streamfunction(x, p, Γ; aggregate=sum)

    nvortices = length(Γ)
    # make sure we have exactly the amount of vortices we need
    @assert length(p) == 2 * nvortices


    distancetovortices = reshape(p, 2, nvortices) .- x
    streampervortex = -(1/2/π) * Γ' .*  
        log.( 
            norm.(eachcol(distancetovortices))  # columnwise length
            )  
    Ψ = aggregate(streampervortex)
    @show Ψ
    Ψ
    
end

export streamfunction


end
