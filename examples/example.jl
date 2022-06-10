# simple demo of vortex dance

using .VortexDynamics2D
using CairoMakie

# set up the initial configuration
N = 5
p = [Vec2D(c) for c in eachcol(rand(2,N)*5)] # initial config
Γ = randn(N) + 2 .* sign.(randn(N)) # circulations

println("Configured the vortex configuration")

# simulate
sol, diffeq = vortexdance( p, Γ, 10.)
println("Completed solving for trajectories")

t = 0:0.01:10
sol_traces = sol(t)


# visualize locations of vortices

	fig = Figure()

	ax = Axis(fig[1,1], 
            xlabel = "x", 
            ylabel = "y", 
            title = "Plot of $N vortices")

	range_of_Γ = Tuple([-1,1] * maximum(abs.(Γ)))
    
    # for each vortex
	for k in eachindex(Γ)
		xs = first.([snapshot[k] for snapshot in sol_traces])
		ys = last.([snapshot[k] for snapshot in sol_traces])

		Makie.lines!(ax,xs,ys,
			legend_label = "Vortex $k",
			color=repeat( Γ[k], inner=(length(t),) ), 
			colormap=:balance,
			colorrange=range_of_Γ
		)
	end
	Colorbar(fig[1,2],limits=range_of_Γ,colormap=:balance)
	
	fig
