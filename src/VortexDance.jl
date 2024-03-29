module VortexDance

using LinearAlgebra:norm
using DifferentialEquations
using StaticArrays

export vortexdance
export biotsavart
export Vec2D
export vortex_ODE

const Vec2D = StaticArrays.SVector{2}


function vortex_ODE(u,p,t)
    out = biotsavart(u, u, p; fieldtype=:velocity)        
    return convert(typeof(u),out)
end

"""
This version of the function creates a return function that can be interpolated
"""
function vortexdance( p_initial, Γ, tspan::Tuple{Real, Real} )
    
    
    diffeq = ODEProblem(
    vortex_ODE,
    p_initial,
    tspan,
    Γ)
    sol = solve( diffeq, Tsit5() )
    
    return sol, diffeq
end

function vortexdance( p_initial, Γ, Tend::Real )

    return vortexdance(p_initial, Γ, (0, Tend) )

end

"""
This version of the function creates an interpolated array
"""
function vortexdance( p_initial, Γ, Trange::AbstractVector )

    # evaluate solution for the endpoints of the range
    sol, diffeq = vortexdance(p_initial, Γ, (first(Trange), last(Trange)) )

    # interpolate solution
    sol_traces = sol(Trange)

    # create a matrix of the solution
    sol_matrix = mapreduce( permutedims, vcat, sol_traces.u )

    return sol_matrix, sol, diffeq

end


function biotsavart( pe::Vec2D, pv::Vector{T}, Γ::Vector; args...) where {T <: Vec2D}
    
    first(biotsavart( [pe,], pv, Γ; args...))
    
end

"""
distances, normals = pdist_normal(p1::AbstractVector{T}, p2::AbstractVector{T}) where {T <: Vec2D}

Compute lengths of all pairs of p1-p2 vectors,
and the normals.

"""
function pdist_normal(p1::AbstractVector{T}, p2::AbstractVector{T}) where {T <: Vec2D}

    # compute normals to pairwise difference vectors between Vec2Ds in pe and pv
    # rotation by 90deg
    S = @SMatrix [0 1; -1 0]

    # compute pairwise differences and rotate them by 90deg
    out = p1 .- permutedims(p2) # pairwise differences
    out = broadcast!( x-> S*x, out, out) # rotate all elements by 90 deg (perps)
    
    # distances between evaluation Vec2D and vortices
    distances = norm.(out)

    return distances, out


end

function biotsavart( pe::AbstractVector{T}, pv::Vector{T}, Γ::Vector; core=1e-12, fieldtype::Union{Symbol,Vector{Symbol}} =[:velocity, :stream]
    ) where {T <: Vec2D}

    # ensure that number of vertices is the same as number of circulations
    @assert length(pv) == length(Γ)

    # ensure fieldtype is either velocity function or stream function
    @assert fieldtype in [:velocity, :stream,:vorticity]

    # ensure that the core is a non-negative number
    @assert core >= 0
    
    distances, normals = pdist_normal(pe, pv)


    if fieldtype == :velocity
    # individual contribution to velocity fieldtype and stream function
    vfs = -(1/2/π) * (Γ' .* normals ) ./ (distances .^ 2 .+ core .^2 )

    vf = similar(pe, VortexDance.Vec2D)
    vf .= sum(vfs, dims=2)

    @assert size(vf) == size(pe)

    return vf
    
    elseif fieldtype==:stream
    streams = (-1/2/π) * Γ' .* log.(distances .+ core)
    
    # sum all contributions
    stream = similar(pe, typeof(streams[1]))
    stream .= sum(streams, dims=2)

    @assert size(stream) == size(pe)
    return stream
elseif fieldtype==:vorticity
    ωs = (1/2/π) * Γ' .* exp.( -0.5.*distances.^2/core.^2 )
    
    # sum all contributions
    ω = similar(pe, typeof(ωs[1]))
    ω .= sum(ωs, dims=2)

    @assert size(ω) == size(pe)
    return ω    
    
    end

    # note to future self: 
    #   once we're trying to implement "max" or "min" 
    #   instead "sum" use findmax/findmin
    
    
end

end

