#region
N = 300;

c1 = [-10.  -10.] .- 2*randn(N,2);
c2 = [ 10.  -10.] .- 2*randn(N,2);
c3 = [-10.   10.] .- 2*randn(N,2);
c4 = [ 10.   10.] .- 2*randn(N,2);

c = [c1;c2;c3;c4];

c = randn(N,2) * 10.0



#Γ = 10.0 * sign.(c[:,2]).*sign.(c[:,1])/N  + randn(4*N,1)/N;
#Γ = ( 10.0 .+ randn(4*N,) ) .* sign.(randn(4*N,))/N
Γ = sign.(randn(N,))*20 + randn(N,)




u0 = c[:]
Γ = Γ[:]/( length(c) ÷ 2)
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
Nv = length(Γ)
pp = @views plot( 
    transpose( sol[1:Nv,:] ), 
    transpose( sol[(Nv+1):end,:]) 
    )
#%xlims!((-15., 15.))
#ylims!((-15., 15.))
@show pp

# using Test
# @testset "VortexDynamics2D.jl" begin
#     # Write your tests here.
# end
