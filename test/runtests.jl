#region
N = 10;
c = randn(N, 2) * 10.0

Γ = sign.(randn(N,)) * 20 + randn(N,)
u0 = Float64.(c[:])
Γ = Float64.(Γ[:]) / N
T = 500.0

sol = VortexDynamics2D.vortexdance(u0,Γ,T)

# using Test
# @testset "VortexDynamics2D.jl" begin
#     # Write your tests here.
# end
