# SargassumBOMB.jl

SargassumBOMB is a pure [Julia](https://julialang.org/) package for computing trajectories, distributions 
and other statistics of pelagic [*Sargassum*](https://oceanexplorer.noaa.gov/facts/sargassum.html). *Sargassum* 
is a brown seaweed which plays a crucial role in the ecosystem of the Sargasso Sea and surrounding areas of the North Atlantic. Since approximately 2011, islands in the Caribbean Ses, beaches in South Florida as well as certain 
regions of western Africa and northern Brazil have been inundated with abnormally large quantities of *Sargassum*. 
The main goal of this package is to provide a performant, extensible and fully open-source toolbox to study *Sargassum* motion.

## Contents

```@contents
Pages = ["index.md"]
```

# Installation 

## Installing Julia

Currently, SargassumBOMB is solely distributed through [GitHub](https://github.com/70Gage70/SargassumBOMB.jl) as a Julia package. It is recommended to install Julia via the [juliaup](https://github.com/JuliaLang/juliaup) version manager.

## Installing SargassumBOMB

In the Julia REPL, run the following

```julia
import Pkg
Pkg.add(url = "https://github.com/70Gage70/SargassumFromAFAI.jl.git")
```

# Quickstart

Want to get started right away? Head to the [first steps](first-steps.md) page of the documentation.


# Package Overview

At the highest level, floating *Sargassum* is modeled as a series of atomic "clumps" connected by non-Hookian springs into a network called a "raft." Each clump obeys a [Maxey-Riley equation](theory); an equation describing the evolution of a small, spherical inertial particle in a background flow. Clumps behave as they would in reality: they beach when reaching land, fracture from one large raft into multiple smaller rafts (and vice cera), and grow and die based on factors like the temperature and nutrient content of the water. Raft trajectories are integrated using Julia's [differential equations ecosystem](https://github.com/SciML/DifferentialEquations.jl), in particular making heavy use of the [callback](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/) functionality.

A significant number of [interpolants](interpolants) are required; these are constructed from raw datasets measuring various parameters in the ocean and atmosphere. SargassumBOMB comes with data covering the year 2018 and adding custom interpolants for other quantities and times is possible. 

SargassumBOMB comes with plotting functionality for trajectories and distributions provided by Julia's [Makie](https://docs.makie.org/stable/) ecosystem.


