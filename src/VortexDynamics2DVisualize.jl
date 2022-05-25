module VortexDynamics2DVisualize

export plotvf

include("VortexDynamics2D.jl")
using .VortexDynamics2D
using CairoMakie: Makie
using SplitApplyCombine:invert
using LinearAlgebra:norm

function plotvf( pv, Γ; N = 21,
    xs = LinRange(-20, 20, N), 
    ys = LinRange(-20, 20, N),
    ax = Makie.current_axis() )
    
    # compute the velocity field
    F(x,y) = VortexDynamics2D.biotsavart([x y],pv,Γ;core=1e-3)
    uv = [ F(x,y) for x in xs, y in ys ] # list of value pairs (pointwise vectors)
    
    speed = norm.(uv)
    
    
    uv = invert(uv) # arrange as a pair of rectangular fields
    
    # clear the axis
    Makie.empty!(ax)
    
    
    
    # plot velocity field colored by the speed
    # logarithms are there to avoid having all visual change happen only very near the vortices themselves
    arrows_handle = Makie.arrows!(xs, ys, 
    uv[1], uv[2], 
    arrowsize = 5, 
    lengthscale = 0.5*sqrt( (xs[2]-xs[1]).^2 + (ys[2]-ys[1]).^2 ), 
    normalize=true,
    colormap = :binary )
    
    # plot point vortices colored by circulation Γ
    scatter_handle = Makie.scatter!(ax, 
    pv[1:(length(pv) ÷ 2)], pv[length(pv)÷2+1:end], 
    color=Γ, colorrange=Tuple(maximum(abs.(Γ)).*[-1,1]),
    colormap=Makie.Reverse(:RdBu_5) )
    
    Makie.limits!(xs[1],xs[end],ys[1],ys[end])
    
    return scatter_handle, arrows_handle
    
    
end

end