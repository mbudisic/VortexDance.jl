using .VortexDynamics2D
using .VortexDynamics2DVisualize

#region
N = 20

c = randn(N, 2) * 10.0
Γ = sign.(randn(N,)) * 20 + randn(N,)
u0 = Float64.(c[:])
Γ = Float64.(Γ[:] ./ N)

sol = vortexdance( u0, Γ, 500)
#endregion

#region
using CairoMakie
using Printf



# animation settings
nframes = 50
time_iterator = range(0, 500, length=nframes)
resolution = 61

t = 0
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1,1])
sc_h, ar_h = plotvf( sol(t),Γ; N=resolution, ax=ax )
Colorbar(fig[1, 2], sc_h, label = "Circulation Γ")
#Colorbar(fig[2, 1], ar_h, label = "Speed")

ax.title = "t = $(@sprintf("%0.3f",t))"
fps=60
for (idx,t) in enumerate(time_iterator)
            ax.title = ("t = $(@sprintf("%0.3f",t))")
            ax
            plotvf( sol(t),Γ; N=resolution, ax=ax )

            println("t = $(@sprintf("%0.3f",t)) done")

            save("vortex_$(@sprintf("%03d",idx)).png",fig)
            #sleep(1/fps)

            
end





#endregion