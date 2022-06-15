# VortexDance.jl [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mbudisic.github.io/VortexDance.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mbudisic.github.io/VortexDance.jl/dev) [![Build Status](https://github.com/mbudisic/VortexDance.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/mbudisic/VortexDance.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![Build Status](https://travis-ci.com/mbudisic/VortexDance.jl.svg?branch=master)](https://travis-ci.com/mbudisic/VortexDance.jl)

Simulate the motion of singular point vortices in a plane. 

Vortices are modeled by positions $p_n \in \mathbb{R}^2$ and circulations $\Gamma_n$, $n=1,2,\dots,N$. Their motion is governed by the Biot--Savart law.
This package defines the corresponding velocity field and stream function, with accessory functions needed for visualizing the outputs.

The motion of vortices is simulated using ![`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl) package.

After cloning the repository, run the Pluto notebook in `examples/pluto_example_motion.jl`. Navigate to the local folder 
```
julia> using Pluto
julia> Pluto.run(notebook="path/to/repo/examples/pluto_example_motion.jl")
```

After playing with the number of vortices and simulation length, you should be able to obtain something like:

![five_vortices](https://user-images.githubusercontent.com/748221/173427025-1abdca3e-d21f-429e-a1c6-846e74b724a7.png)
