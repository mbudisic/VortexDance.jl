module VortexDance

using LinearAlgebra:norm
using DifferentialEquations
using StaticArrays

export vortexdance
export biotsavart
export Vec2D

const Vec2D{T} = StaticArrays.SVector{2,T}

"""
This version of the function creates a return function that can be interpolated
"""
function vortexdance( p_initial, Γ, tspan::Tuple{Real, Real} )
    
    velocityfield(u, p, t) = biotsavart(u, u, p; field=:velocity)
        
    diffeq = ODEProblem(
    velocityfield,
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
function vortexdance( p_initial, Γ, T::AbstractVector )

    # evaluate solution for the endpoints of the range
    sol, diffeq = vortexdance(p_initial, Γ, (first(T), last(T)) )

    # interpolate solution
    sol_traces = sol(T)

    # create a matrix of the solution
    sol_matrix = mapreduce( permutedims, vcat, sol_traces.u )

    return sol_matrix, sol, diffeq

end


function biotsavart( pe::Vec2D, pv::Vector{T}, Γ::Vector; args...) where {T <: Vec2D}
    
    retvalues =  biotsavart( [pe], pv, Γ; args...)
    return first.( (retvalues...) )

end

"""


"""
function biotsavart( pe::AbstractVector{T}, pv::Vector{T}, Γ::Vector; core=1e-12, field=(:velocity, :stream)
    ) where {T <: Vec2D}
    
    # ensure that number of vertices is the same as number of circulations
    @assert length(pv) == length(Γ)

    # ensure that the core is a non-negative number
    @assert core >= 0
    
    # compute normals to pairwise difference vectors between Vec2Ds in pe and pv
    S = @SMatrix [0 1; -1 0]
    normals = broadcast( 
        x->S*x, 
        reshape( view( pe, :, : ), : , 1 ) .- reshape( view( pv, :, :), 1, :  ) 
        )
    
    # distances between evaluation Vec2D and vortices
    distances = norm.(normals)
    

    retfields = ()

    if :velocity in field
    # individual contribution to velocity field and stream function
    vfs = -(1/2/π) * (Γ' .* normals ) ./ (distances .^ 2 .+ core .^2 )
    vf = vec(  sum(vfs, dims=2) )

    retfields = (retfields..., vf)
    end

    if :stream in field
    streams = (-1/2/π) * Γ' .* log.(distances .+ core)
    
    # sum all contributions
    stream = vec( sum(streams, dims=2) )

    retfields = (retfields..., stream)
    end

    # note to future self: 
    #   once we're trying to implement "max" or "min" 
    #   instead "sum" use findmax/findmin
    
    return retfields
    
end

end

