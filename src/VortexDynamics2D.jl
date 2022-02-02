module VortexDynamics2D

"""
    streamfunction(x; p, Γ)

Evaluate the stream function of a collection of 2d vortices with positions
in vector `p` and circulations in vector `Γ` using the Biot--Savart law.


The length of vector `Γ` determines the number of vortices `N`.
The length of vector `p` should be exactly `2N` with coordinates of vortices 
interlaced, i.e., `x,y,x,y,...` for 
"""
function streamfunction(x, p, Γ)

    
    2x + 3y
    
end

export streamfunction


end
