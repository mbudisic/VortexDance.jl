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
	Pkg.add("PlutoUI")
	using PlutoUI	
md"""
Cell hidden - packages here are just needed for the notebook.
"""	
end

# ╔═╡ 92a981be-1161-4cc8-a404-be0bdeab2903
begin
	Pkg.add(
		url="git@github.com:mbudisic/VortexDance.jl.git",
		rev="master"
	) # the github url can be replaced by a local path to the repository
	using VortexDance
end

# ╔═╡ 0d8aded3-5ec4-49f8-8eb4-f49cdde86624
md"""
# Simple demo of vortex dance

This is the demo of the `VortexDance.jl` package's basic functionality --- solve for the motion of N vortices in the plane evolving under the Biot--Savart law.

This is a Pluto.jl notebook, which means that the cells do not have to be in the order of execution by `julia`. To see the "correct" order of cells, open this file in an editor.
"""

# ╔═╡ a766f8e2-e011-4689-b5ef-ff1387c018f4
md"""
First, let' see what the final product works.

The slider defines how many vortices are there, and their positions and circulations are chosen randomly.

Vortices: $N=$ $(@bind N PlutoUI.Slider(1:20,show_value=true))

Duration: $T=$ $(@bind T PlutoUI.Slider(10:5:100,show_value=true))

Note: if the vortices are randomly placed very close to each other, the ODE becomes more "stiff", and the solver takes longer to produce the solution.
"""

# ╔═╡ c1c57edd-b02b-4891-8b7c-f87df83c4afe
md"""
### Going through the code
"""

# ╔═╡ d8275ac4-f301-45f4-92f0-791e76f14c23
md"""
Here we load the `VortexDance.jl` package from the `git` repository. This syntax allows us *not* to have to register the package with the Julia package system (since we're still developing it) AND allows us to pull a particular branch of the `git` repository.	
	"""

# ╔═╡ bbfabf36-9659-4790-80bd-2b980248cf82
md"""
The `VortexDance` package uses a custom `Vec2D` datatype for locations of vortices. `Vec2D` is simply a `StaticArrays` vector of two elements, which allows for a more straightforward codes downstream.

We initialize a vector of vortices by generating a random matrix, and then converting each column into a `Vec2D`.

(Remember: `N` was defined by the slider at the beginning of the notebook.)
"""

# ╔═╡ b2f950bf-6966-4f04-8399-1fa5ccd2c112
x = [Vec2D(c) for c in eachcol(rand(2,N)*8)]; # semicolon prevents output

# ╔═╡ 94f17601-cdec-4a08-84b5-916d386978fa
md"""
Circulations $\Gamma$ are chosen as fluctuations around $\pm 2$.
"""

# ╔═╡ bda16484-a8d0-4c34-be10-a37870444017
Γ = randn(N) + 2 .* sign.(randn(N))

# ╔═╡ 674ff49a-a8e5-4a0c-b4f1-6cdb25324b4b
md"""
Function `VortexDance.vortexdance(...)` takes in initial locations, circulations, duration, and produces the trajectories. Remember: `T` was defined by the slider at the beginning of the notebook.
"""

# ╔═╡ e31220b4-5d2d-4006-b302-773ce8bdf6d9
sol_matrix, sol, diffeq  = vortexdance( x, Γ, 0:0.01:T);

# ╔═╡ 6feabdee-f7cc-4fda-af66-de4962fcc38a
md"""
The first output is the solution interpolated at the range of values passed in. It's produced only if the time is specified as a vector of values. This is in the "matrix of points" format - rows are timesteps, columns correspond to vortices.

Second output is a solution structure that can be evaluated at any point in time. e.g. running `sol(π)=`
"""

# ╔═╡ 0153f562-0e64-4906-bd47-42e56aa912a1
display(sol(π))

# ╔═╡ 0d3afbe2-d998-4033-a486-a40eccd65094
md"""
Third output is the differential equation structure, that stores all kinds of info about a problem, and could be re-solved using a different numerical solver.
"""

# ╔═╡ 86eda4c7-d596-4f1c-8105-68930c94ac07
md"""
The next block of code produces the figure of the vortex solution. Each vortex trace is plotted by its own `Makie.lines` command. 
"""

# ╔═╡ 4ae17fc3-1779-454d-9061-2575e7adfca5
# invoking vortex_figure in a cell displays the output of this block
vortex_figure = with_theme(theme_minimal()) do
	
	fig = Figure(resolution=(400,400))

	# figure has two Axes: "main" ....
	ax = Axis(fig[1,1], 
            xlabel = "x", 
            ylabel = "y", 
            title = "$N vortices, t ∈ [0,$T]",
			aspect=1)

	# ... and "colorbar"
	range_of_Γ = Tuple([-1,1] * maximum(abs.(Γ)))
    Colorbar(fig[1,2],
		limits=range_of_Γ,
		colormap=:balance,
		label="Circulation Γ")

    # iterate over vortex trajectories
	# this syntax assigns variable name to index and value of each element
	for (k, vortex) in enumerate(eachcol(sol_matrix))

		# extract first/last element (=x,y coordinates) for each vortex trajectory
		xs = first.(sol_matrix[:,k])
		ys = last.(sol_matrix[:,k])

		# make the color vector of the same size
		mycolor = fill( Γ[k], size(sol_matrix,1) )

		# display a line
		Makie.lines!(ax,xs, ys,
			legend_label = "Vortex $k",
			color=mycolor, 
			colormap=:balance,
			colorrange=range_of_Γ
		)
	end

	return fig
end;

# ╔═╡ 0a22777a-fbaa-4f15-be84-95fbef016050
md"""
The following code simply outputs the produced figure to the browser. Later, we're going to see how this `fig` is generated.

$(vortex_figure)

Vortex color is determined by their circulation.

"""

# ╔═╡ fad93175-d66f-4d08-86ab-e62311b657a6


# ╔═╡ Cell order:
# ╟─0d8aded3-5ec4-49f8-8eb4-f49cdde86624
# ╟─a766f8e2-e011-4689-b5ef-ff1387c018f4
# ╟─0a22777a-fbaa-4f15-be84-95fbef016050
# ╟─c1c57edd-b02b-4891-8b7c-f87df83c4afe
# ╠═20b1ef12-e8f0-11ec-0518-ef96740735a8
# ╟─d8275ac4-f301-45f4-92f0-791e76f14c23
# ╠═92a981be-1161-4cc8-a404-be0bdeab2903
# ╟─bbfabf36-9659-4790-80bd-2b980248cf82
# ╠═b2f950bf-6966-4f04-8399-1fa5ccd2c112
# ╟─94f17601-cdec-4a08-84b5-916d386978fa
# ╠═bda16484-a8d0-4c34-be10-a37870444017
# ╟─674ff49a-a8e5-4a0c-b4f1-6cdb25324b4b
# ╠═e31220b4-5d2d-4006-b302-773ce8bdf6d9
# ╟─6feabdee-f7cc-4fda-af66-de4962fcc38a
# ╟─0153f562-0e64-4906-bd47-42e56aa912a1
# ╟─0d3afbe2-d998-4033-a486-a40eccd65094
# ╟─86eda4c7-d596-4f1c-8105-68930c94ac07
# ╠═4ae17fc3-1779-454d-9061-2575e7adfca5
# ╠═fad93175-d66f-4d08-86ab-e62311b657a6
