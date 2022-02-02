module VortexDynamics2D

using LinearAlgebra
using Infiltrator

"""
    streamfunction(x; p, Γ)

Evaluate the stream function of a collection of 2d vortices with positions
in vector `p` and circulations in vector `Γ` using the Biot--Savart law.


The length of vector `Γ` determines the number of vortices `N`.
The length of vector `p` should be exactly `2N` with coordinates of vortices 
interlaced, i.e., `x,y,x,y,...` for 
"""
function streamfunction(x, p, Γ)

    # make sure we have exactly the amount of vortices we need
    
    nvortices = length(Γ)
    
    @assert length(p) == 2 * nvortices

    dmat = reshape(p, 2, nvortices) .- x
    streampervortex = - Γ' .* log.( norm.(eachcol(dmat)) ) /2 /π

    Ψ = maximum(streampervortex)

    Ψ
    
end

export streamfunction


end
