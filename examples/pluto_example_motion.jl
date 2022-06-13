### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 20b1ef12-e8f0-11ec-0518-ef96740735a8
# simple demo of vortex dance
begin
	using Revise
	using Pkg
	using CairoMakie
end


# ╔═╡ 92a981be-1161-4cc8-a404-be0bdeab2903
begin
	Pkg.add(
		url="git@github.com:mbudisic/VortexDynamics2D.jl.git",
		rev="TwoColumnState"
	)
	using VortexDynamics2D
end

# ╔═╡ bafe96b0-2207-4533-ba11-09de45d2e0d5
begin
	Pkg.add("PlutoUI")
	using PlutoUI
end

# ╔═╡ 0d8aded3-5ec4-49f8-8eb4-f49cdde86624
md"""
# Simple demo of vortex dance
"""

# ╔═╡ bbfabf36-9659-4790-80bd-2b980248cf82
@bind N PlutoUI.Slider(3:20,show_value=true)

# ╔═╡ b2f950bf-6966-4f04-8399-1fa5ccd2c112

begin
# set up the initial configuration
x = [Vec2D(c) for c in eachcol(rand(2,N)*5)] # initial config
Γ = randn(N) + 2 .* sign.(randn(N)) # circulations
end


# ╔═╡ 4d5bcded-a614-4874-9e8b-ad17e55a19a1

# simulate
begin
	sol, diffeq = vortexdance( x, Γ, 10.)

	t = 0:0.01:10
	sol_traces = sol(t)	
end


# ╔═╡ 4ae17fc3-1779-454d-9061-2575e7adfca5
# visualize locations of vortices
begin
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
			color=repeat( [Γ[k]], inner=(length(t),) ), 
			colormap=:balance,
			colorrange=range_of_Γ
		)
	end
	Colorbar(fig[1,2],limits=range_of_Γ,colormap=:balance)
	
	fig
end

# ╔═╡ 63615a2e-57f2-4f6d-968a-57bf28e8ce4a


# ╔═╡ Cell order:
# ╟─0d8aded3-5ec4-49f8-8eb4-f49cdde86624
# ╠═20b1ef12-e8f0-11ec-0518-ef96740735a8
# ╠═92a981be-1161-4cc8-a404-be0bdeab2903
# ╠═bafe96b0-2207-4533-ba11-09de45d2e0d5
# ╠═bbfabf36-9659-4790-80bd-2b980248cf82
# ╠═b2f950bf-6966-4f04-8399-1fa5ccd2c112
# ╠═4d5bcded-a614-4874-9e8b-ad17e55a19a1
# ╠═4ae17fc3-1779-454d-9061-2575e7adfca5
# ╠═63615a2e-57f2-4f6d-968a-57bf28e8ce4a
