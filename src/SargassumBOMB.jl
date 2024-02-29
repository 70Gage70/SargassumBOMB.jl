module SargassumBOMB

# core functionality
using OrdinaryDiffEq, SargassumFromAFAI, NearestNeighbors, Interpolations
using Unitful, Dates
using LinearAlgebra: norm, ⋅

# i/o
using MAT, NetCDF, JLD2

# probability/statistics
using StatsBase, Distributions 
using Random: seed!

# plotting/geography
using Makie, CairoMakie, GeoMakie
using SargassumColors
using GeoMakie, GeoMakie.GeoJSON, GeoDatasets, GeoInterface

# printing
using Crayons.Box
using Latexify

# optimization 
using Metaheuristics
using QuasiMonteCarlo

###############################################################

include("utils.jl")
export vec2range

include("coordinates.jl")
export UNITS, EARTH_RADIUS, EquirectangularReference, EQR, sph2xy, xy2sph

include("time.jl")
export T_REF, datetime2time, ymw2time, ymwspan2weekspan

include(joinpath(@__DIR__, "..", "interpolants", "itps-core.jl"))
export GriddedField, InterpolatedField
export add_spatial_dimension!, add_temporal_dimension!, add_field!, ranges_increasing!, sph2xy!, add_derivatives!

include(joinpath(@__DIR__, "..", "interpolants", "itps-construct.jl"))
export WATER_ITP, WIND_ITP, WAVES_ITP, STOKES_ITP, LAND_ITP, TEMPERATURE_ITP, NUTRIENTS_ITP
export itps_default_construct, itps_default_assign

include(joinpath(@__DIR__, "biology.jl"))
export AbstractGrowthDeathModel, ImmortalModel, BrooksModelParameters, brooks_dSdt_clump, brooks_dSdt_raft, BrooksModel

include("geography.jl")
export AbstractLand, NoLand, Land

include("raft-parameters.jl")
export ClumpParameters
export SpringParameters, ΔL, spring_force, BOMB_k
export InitialConditions
export AbstractConnections, ConnectionsNone, ConnectionsFull, ConnectionsRadius, ConnectionsNearest, form_connections!
export RaftParameters

include("physics.jl")
export Raft!, Leeway!

include("control.jl")
export n_clumps, clump_i, com, cb_update, cb_land, cb_growth_death, cb_connections, kill!, grow!

include("trajectories.jl")
export Trajectory, time_slice, RaftTrajectory, uniformize, bins

include("main.jl")
export simulate

include("optimization.jl")
export OPTIMIZATION_PARAMETER_NAMES, LossFunction, OptimizationParameter, BOMBOptimizationProblem, loss, optimize!, sample!, save_bop, load_bop

include(joinpath(@__DIR__, "..", "plotting", "plotting-core.jl"))
export default_fig, geo_axis, land!, data_legend!, trajectory!, trajectory_hist!, plot

include(joinpath(@__DIR__, "..", "plotting", "plotting-itp.jl"))
export check_land, check_itp

export length, show, iterate # various Base extensions

# initialize interpolants
function __init__()
    itps_default_construct()
    itps_default_assign()
end

import PrecompileTools

PrecompileTools.@compile_workload begin
    itps_default_construct()
    itps_default_assign()

    tspan = (0.0, 5.0)
    ics = InitialConditions(tspan, range(-55.0, -50.0, length = 5), range(5.0, 10.0, length = 5))
    clumps = ClumpParameters()
    springs = SpringParameters(k -> 0.1, 100.0)
    connections = ConnectionsNearest(10)
    gd_model = BrooksModel()
    land = Land()
    

    rp = RaftParameters(
        ics = ics,
        clumps = clumps,
        springs = springs,
        connections = connections,
        gd_model = gd_model,
        land = land
    )

    sol = simulate(rp, showprogress = false)
end

end # module

