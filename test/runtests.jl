#region
N = 4;

c1 = [-10.  -10.] .- 2*randn(N,2);
c2 = [ 10.  -10.] .- 2*randn(N,2);
c3 = [-10.   10.] .- 2*randn(N,2);
c4 = [ 10.   10.] .- 2*randn(N,2);

c = [c1;c2;c3;c4];



Γ = 10.0 * sign.(c[:,2]).*sign.(c[:,1])/N  + randn(4*N,1)/N;

Γ = ( 10.0 .+ randn(4*N,) ) .* sign.(randn(4*N,))/N





u0 = c[:]
Γ = Γ[:]
#endregion

∇Ψ = VortexDynamics2D.biotsavart(c,c,Γ)

vf(u,p,t) = VortexDynamics2D.biotsavart(u,u,p)

using DifferentialEquations

tspan = (0, 500.0)

prob = ODEProblem(
    vf,
    u0,
    tspan,
    Γ)

sol = solve(prob,Tsit5(),
    abstol=1e-6, reltol=1e-6,saveat=0.01)

using Plots
M = length(Γ)
pp = plot( sol[1,:], sol[1+1,:])
for k=3:2:M
    plot!(pp, sol[k,:], sol[k+1,:])
end
#%xlims!((-15., 15.))
#ylims!((-15., 15.))
@show pp

# using Test
# @testset "VortexDynamics2D.jl" begin
#     # Write your tests here.
# end
