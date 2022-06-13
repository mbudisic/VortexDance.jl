using RectiGrids
using StaticArrays

XYgrid = RectiGrids.grid(SVector{2,}, -20.:20, -40.:40)

vf(xy) = VortexDance.biotsavart(xy, u0, Γ; core=1e-3)

VF = vf.(XYgrid)

XY = invert(XY)
VF = invert(VF)

scale = 1000
Plots.quiver( XY[1][:], XY[2][:], quiver=(scale.*VF[1][:], scale.*VF[2][:]))