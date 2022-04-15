module VortexDynamics2D

using LinearAlgebra
using StaticArrays

"""
    biotsavart(pv, pe, Γ; aggregate=sum)

Evaluate the stream function of a collection of 2d vortices with positions
in vector `p` and circulations in vector `Γ` using the Biot--Savart law.

`x` is the position of evaluation
The length of vector `Γ` determines the number of vortices `N`.
The length of vector `p` should be exactly `2N` with coordinates of vortices 
stacked, i.e., `x1,...,xN,y1,...,yN` 

setting the keyword argument `aggregate` to something that is not `sum`
uses a different aggregation function across vortices. e.g. `aggregate=maximum` 
will use the strongest contribution, not the sum of contributions. This is nonphysical 
but can be of interest.
"""
function biotsavart(pe, pv, Γ; aggregate=:sum, core=1.0e-6)

    #TODO Convert to non-allocating
    #TODO Add multiple-dispatch depending on the shape of inputs pv pe

    Nv = length(pv) ÷ 2
    Ne = length(pe) ÷ 2

    pvx = @views pv[1:Nv]
    pvy = @views pv[(Nv+1):end]

    pex = @views pe[1:Ne]
    pey = @views pe[(Ne+1):end]

    # ensure that number of vertices is the same as number of circulations
    @assert size(pvx) == size(Γ)

    D = zeros(Ne, Nv, 2)
    Lsq = zeros(Ne, Nv, 1)

    vectors_normal_to_pairs!(D, Lsq, [pex pey], [pvx pvy], core)
    # compute the Biot--Savart
    perp_gradΨ = -(transpose(Γ) / 2 / π) .* D ./ Lsq

    # add all vortex influences  size(V) = #evaluation points x 1 x 2
    V = sum(perp_gradΨ, dims=2)
    V = V[:]


end

"""
Computes perps to p1-p2 and the euclidean-squared norm.


both p1 and p2 are assumed to be N1x2, N2x2 matrices
D is N1 x N2 x 2
L is N1 x N2 x 1

If optional argument is passed, distances shorter than core_dist are 
set to core_dist, and the associated perp vectors are set to zero vectors
"""
function vectors_normal_to_pairs!(D, Lsq, p1, p2, core_dist=0)
    @views D[:, :, 1] .= p1[:, 2] .- transpose(p2[:, 2])
    @views D[:, :, 2] .= -(p1[:, 1] .- transpose(p2[:, 1]))

    sum!(Lsq, D .^ 2) # compute euclidean norm for each D(i,j,:)

    sel = @. Lsq <= core_dist^2
    @. D[sel, :] = 0
    @. Lsq[sel] = core_dist^2

end

end