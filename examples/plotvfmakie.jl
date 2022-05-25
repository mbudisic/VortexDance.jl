using CairoMakie

hfig = Figure(resolution = (800, 800))
ax =Axis(f[1, 1], backgroundcolor = "black")
N = 61
xs = LinRange(-20, 20, N)
ys = LinRange(-20, 20, N)

F(x,y) = VortexDynamics2D.biotsavart([x y],u0,Γ;core=1e-3)

uv = [ F(x,y) for x in xs, y in ys ]
strength = norm.(uv)
uv = invert(uv)

arrows(xs, ys, log.(strength).*uv[1]./strength, log.(strength).*uv[2]./strength, arrowsize = 5, 
    lengthscale = 0.25, arrowcolor=strength,
    axis = (limits = (xs[1],xs[end],ys[1],ys[end]),) )#, arrowcolor = strength, linecolor = strength)
