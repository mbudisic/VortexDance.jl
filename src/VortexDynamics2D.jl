module VortexDynamics2D

using LinearAlgebra:norm
using DifferentialEquations
using StaticArrays

export vortexdance
export biotsavart
export vectornormal

const Point{T} = StaticArrays.SVector{2,T}

"""



"""
function vortexdance( p_initial, Γ, T )

    vf(x) = biotsavart(x,x,Γ)


    velocityfield(u, p, t) = VortexDynamics2D.biotsavart(u, u, p)
    
    tspan = (0, T)

    diffeq = ODEProblem(
    velocityfield,
    p_initial,
    tspan,
    Γ)
    sol = solve( diffeq, AutoTsit5(Rosenbrock23()) )

    return sol, diffeq
end

function biotsavart( pe::Vector{T}, pv::Vector{T}, Γ::Vector; core=1e-12) where {T <: Point}

	f(x) = biotsavart(x, pv, Γ; core)
	return f.(pe)

end

"""


"""
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

