#region
N = 20

c = randn(N, 2) * 10.0
Γ = sign.(randn(N,)) * 20 + randn(N,)
u0 = Float64.(c[:])
Γ = Float64.(Γ[:] ./ N)

sol = VortexDynamics2D.vortexdance( u0, Γ, 500)
#endregion

#region
using CairoMakie
using Printf



# animation settings
nframes = 250
framerate = 30
time_iterator = range(0, 500, length=nframes)
resolution = 61

t = 0
fig = Figure()
ax = Axis(fig[1,1])
sc_h, ar_h = VortexDynamics2D.plotvf( sol(t),Γ; N=resolution, ax=ax )
Colorbar(fig[1, 2], sc_h, label = "Circulation Γ")
#Colorbar(fig[2, 1], ar_h, label = "Speed")

ax.title = "t = $(@sprintf("%0.3f",t))"

record(fig, "test.mp4", time_iterator;
        framerate = framerate) do t

            ax.title = ("t = $(@sprintf("%0.3f",t))")
            ax
            VortexDynamics2D.plotvf( sol(t),Γ; N=resolution, ax=ax )

            println("t = $(@sprintf("%0.3f",t)) done")


            
        end





#endregion