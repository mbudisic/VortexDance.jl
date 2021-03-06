# VortexDance.jl [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mbudisic.github.io/VortexDance.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mbudisic.github.io/VortexDance.jl/dev) [![Build Status](https://github.com/mbudisic/VortexDance.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/mbudisic/VortexDance.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![Build Status](https://travis-ci.com/mbudisic/VortexDance.jl.svg?branch=master)](https://travis-ci.com/mbudisic/VortexDance.jl)

Simulate the motion of singular point vortices in a plane. 

Vortices are modeled by positions and circulations $$p_n \in \mathbb{R}^2, \Gamma_n,\quad n=1,2,\dots,N.$$ Their motion is governed by the Biot--Savart law.
This package defines the corresponding velocity field and stream function, with accessory functions needed for visualizing the outputs.

The motion of vortices is simulated using ![`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl) package.

After cloning the repository, run the Pluto notebook in `examples/pluto_example_motion.jl`. Navigate to the local folder 
```
julia> using Pluto
julia> Pluto.run(notebook="path/to/repo/examples/pluto_example_motion.jl")
```

After playing with the number of vortices and simulation length, you should be able to obtain something like:

![vortex-dance-demo](https://user-images.githubusercontent.com/748221/173864068-c8671098-a962-4602-a052-ab6520af4ed2.gif)
