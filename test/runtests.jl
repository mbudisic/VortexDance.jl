N = 2;

c1 = [-3.  -2.] .- rand(N,2);
c2 = [ 3.  -2.] .- rand(N,2);
c3 = [-3.   2.] .- rand(N,2);
c4 = [ 3.   2.] .- rand(N,2);

c = [c1;c2;c3;c4];

Γ = 10.0 * sign.(c[:,2]).*sign.(c[:,1])/N + randn(4*N,1)/N;

c = c[:]
Γ = Γ[:]

∇Ψ = VortexDynamics2D.biotsavart(c,c,Γ)

vf(u,p,t) = VortexDynamics2D.biotsavart(u,u,Γ)

using DifferentialEquations

u0 = c[:]
tspan = (0, 1000.0)

prob = ODEProblem(
    (u,p,t) -> VortexDynamics2D.biotsavart(u,u,Γ),
    u0,
    tspan)

sol = solve(prob,AutoTsit5(Rosenbrock23()),reltol=1e-3,abstol=1e-3,saveat=0.1)

M = length(Γ)
pp = plot( sol[1,:], sol[1+M,:])
for k=2:M
    plot!(pp, sol[k,:], sol[k+M,:])
end
@show pp

# using Test
# @testset "VortexDynamics2D.jl" begin
#     # Write your tests here.
# end
