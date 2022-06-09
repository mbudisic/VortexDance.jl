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

# ╔═╡ 773d1e80-e81c-11ec-17ac-e35e8efb7de2
using StaticArrays

# ╔═╡ 7be5e280-156e-416d-965a-71bbca10c8a5
using LinearAlgebra:norm

# ╔═╡ 40f5344a-051b-422c-bd56-97d27c52bffe
using PlutoUI

# ╔═╡ acf4eb99-7c8e-4c28-9725-98548bfba47b
const Point{T} = StaticArrays.SVector{2,T}

# ╔═╡ de576bd2-603b-479f-945b-e59ea4834161
v = Point([0.1 0.2])

# ╔═╡ ceb5802c-0db2-4e48-b939-8473bbfc955d
begin
	ps = Vector([Point([1, 0]), Point([0.5, 0.5]), Point([-1, 0])])
	Γ = [-1,  1.5,  -0.5]
	typeof(ps)
end

# ╔═╡ 5015418b-1301-4c99-a476-a1ed0c4fe2a6
S = @SMatrix [0 1; -1 0]

# ╔═╡ 86a36ad6-47d2-4ab6-8802-10f7dcb61534
d = ([v] .- ps)

# ╔═╡ 99a9b8c5-2ddf-4b2c-a0fb-52c6dccff164
n = broadcast( x -> S*x, d)

# ╔═╡ 927c9765-dcec-4930-b5e8-0e8a9c2faefa
norms = norm.(n)

# ╔═╡ 9cd135bf-abc5-4e48-9df9-ff4e79aa75ea


# ╔═╡ b3e3a509-6b7e-4117-ac53-44abb79a1d82


# ╔═╡ 0bc563fa-96b2-4d58-ad39-15773204f28c
corebind = @bind corexp Slider(-9:0)

# ╔═╡ 32895524-0fef-40c0-9081-33135bd75b2d
begin
coreval = (10.)^corexp
sel = norms .> coreval
coreval,sel
end

# ╔═╡ a5ffe661-1433-41c1-aeaa-b43ac01c0996
sum( -(1/2/π) * (Γ[sel] .* n[sel]) ./ norms[sel].^2 )

# ╔═╡ bdb28cff-5f41-4457-b542-b72cd1112644
function biotsavart( pe::Point, pv::Vector{T}, Γ::Vector; core=1e-12) where {T <: Point}

	# ensure that number of vertices is the same as number of circulations
    @assert length(pv) == length(Γ)

    # compute normal to distance to point
    S = @SMatrix [0 1; -1 0]
    normals = broadcast( x->S*x, [pe] .- pv )

    # distances between evaluation point and vortices
    distances = norm.(normals)

    # avoid singularity by computing only outside cores
    sel = distances .> core

    # individual contribution to velocity field and stream function
    vfs = -(1/2/π) * (Γ[sel] .* normals[sel]) ./ (distances[sel] .^ 2)
    streams = (-1/2/π) * Γ[sel] .* log.(distances[sel])

    vf = sum(vfs)
    stream = sum(streams)

    return vf, stream





end

# ╔═╡ bd6960ad-cbad-4ea2-a918-ed5b867773da
biotsavart(v, ps, Γ; core=coreval)

# ╔═╡ 82e65e35-591f-4f60-b780-f41fadaa425b
begin
	vf(x) = biotsavart(x, ps, Γ)[1]
	stream(x) = biotsavart(x, ps, Γ)[2]
end

# ╔═╡ 37fce241-18eb-46ee-a59e-f232ef680120
vf.(ps)

# ╔═╡ 380ca39b-750f-4e75-81f5-e9a3c9b24f4b
stream.(ps)

# ╔═╡ d4fca53c-69a7-434b-88a6-5c7c1a3c5d8e
typeof(ps)

# ╔═╡ 4dd4e8fa-6cc1-4299-a865-19a577415617
typeof(v)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
PlutoUI = "~0.7.39"
StaticArrays = "~1.4.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "0f4e115f6f34bbe43c19751c90a38b2f380637b9"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.3"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "383a578bdf6e6721f480e749d503ebc8405a0b22"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.6"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═773d1e80-e81c-11ec-17ac-e35e8efb7de2
# ╠═acf4eb99-7c8e-4c28-9725-98548bfba47b
# ╠═de576bd2-603b-479f-945b-e59ea4834161
# ╠═ceb5802c-0db2-4e48-b939-8473bbfc955d
# ╠═5015418b-1301-4c99-a476-a1ed0c4fe2a6
# ╠═86a36ad6-47d2-4ab6-8802-10f7dcb61534
# ╠═99a9b8c5-2ddf-4b2c-a0fb-52c6dccff164
# ╠═7be5e280-156e-416d-965a-71bbca10c8a5
# ╠═927c9765-dcec-4930-b5e8-0e8a9c2faefa
# ╠═32895524-0fef-40c0-9081-33135bd75b2d
# ╠═9cd135bf-abc5-4e48-9df9-ff4e79aa75ea
# ╠═b3e3a509-6b7e-4117-ac53-44abb79a1d82
# ╠═40f5344a-051b-422c-bd56-97d27c52bffe
# ╠═0bc563fa-96b2-4d58-ad39-15773204f28c
# ╠═a5ffe661-1433-41c1-aeaa-b43ac01c0996
# ╠═bdb28cff-5f41-4457-b542-b72cd1112644
# ╠═bd6960ad-cbad-4ea2-a918-ed5b867773da
# ╠═82e65e35-591f-4f60-b780-f41fadaa425b
# ╠═37fce241-18eb-46ee-a59e-f232ef680120
# ╠═380ca39b-750f-4e75-81f5-e9a3c9b24f4b
# ╠═d4fca53c-69a7-434b-88a6-5c7c1a3c5d8e
# ╠═4dd4e8fa-6cc1-4299-a865-19a577415617
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
