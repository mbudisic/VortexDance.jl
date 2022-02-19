using VortexDynamics2D


N = 4;

c1 = [-8.  -2.] .- rand(N,2);
c2 = [ 8.  -2.] .- rand(N,2);
c3 = [-8.   2.] .- rand(N,2);
c4 = [ 8.   2.] .- rand(N,2);

c = [c1;c2;c3;c4];

Γ = -10.0 * sign.(c[:,2]).*sign.(c[:,1])/N + randn(4*N,1)/N; 

∇Ψ = VortexDynamics2D.biotsavart(c,c,Γ)

vf(u,p,t) = VortexDynamics2D.biotsavart(u,u,Γ[:])

using DifferentialEquations

u0 = c[:]
tspan = (0, 20.0)

prob = ODEProblem(vf, u0, tspan)

sol = solve(prob, saveat=0.1)


# using Test
# @testset "VortexDynamics2D.jl" begin
#     # Write your tests here.
# end
