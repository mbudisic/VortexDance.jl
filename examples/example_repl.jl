includet("../src/VortexDance.jl")

using .VortexDance
using Test
using RecursiveArrayTools
using BenchmarkTools
N = 3; pe = VectorOfArray( [VortexDance.Vec2D(c) for c in eachcol(rand(2,N))] );
M = 4; pv = VectorOfArray( [VortexDance.Vec2D(c) for c in eachcol(rand(2,M))] );
Γ = collect(0:M-1)

a = VortexDance.Vec2D([0.1, 0.2])

biotsavart(pv, pv, Γ; fieldtype = :velocity)
biotsavart(pe, pv, Γ; fieldtype = :velocity)

@test pv isa typeof( biotsavart(pv, pv, Γ; fieldtype = :velocity) )
@test pe isa typeof( biotsavart(pe, pv, Γ; fieldtype = :velocity) )

biotsavart(a, pv,Γ; fieldtype = :velocity)


@btime biotsavart(pv, pv, Γ; fieldtype = :velocity)
@btime biotsavart(pe, pv, Γ; fieldtype = :velocity)
