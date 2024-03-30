module SargassumBOMB

# core functionality
using OrdinaryDiffEq, SargassumFromAFAI, NearestNeighbors, Interpolations
using Unitful, Dates
using LinearAlgebra: norm, ⋅
using ArgCheck

# i/o
using MAT, NetCDF, JLD2, RemoteFiles

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
export n_clumps, clump_i, com, vec2range

include("coordinates.jl")
export UNITS, EARTH_RADIUS, EquirectangularReference, EQR
export sph2xy, xy2sph

include("time.jl")
export T_REF, datetime2time, time2datetime, ymw2time, ymwspan2weekspan

include(joinpath(@__DIR__, "..", "interpolants", "itps-core.jl"))
export GriddedField, InterpolatedField
export add_spatial_dimension!, add_temporal_dimension!, add_field!, ranges_increasing!, sph2xy!, add_derivatives!

include(joinpath(@__DIR__, "..", "interpolants", "itps-construct.jl"))
export WATER_ITP, WIND_ITP, WAVES_ITP, STOKES_ITP, LAND_ITP, TEMPERATURE_ITP, NUTRIENTS_ITP
export itps_default_construct, itps_default_assign

include("land.jl")
export AbstractLand, NoLand, Land

include("ics.jl")
export InitialConditions

include("springs.jl")
export AbstractSpring, HookeSpring, BOMBSpring, ΔL, spring_force
export AbstractConnections, ConnectionsNone, ConnectionsFull, ConnectionsRadius, ConnectionsNearest, form_connections!

include("growth-death.jl")
export AbstractGrowthDeathModel, ImmortalModel, BrooksModelParameters, BrooksModel

include("rafts-clumps.jl")
export ClumpParameters, RaftParameters

include("physics.jl")
export Raft!, Leeway!

include("control.jl")
export cb_update, cb_land, cb_growth_death, cb_connections, kill!, grow!

include("trajectories.jl")
export Trajectory, time_slice, RaftTrajectory, uniformize, bins

include("main.jl")
export simulate

include("optimization.jl")
export OPTIMIZATION_PARAMETER_NAMES, LossFunction, OptimizationParameter, BOMBOptimizationProblem
export optimizable, optimize!, sample!

include(joinpath(@__DIR__, "..", "plotting", "plotting-core.jl"))
export default_fig, geo_axis, land!, data_legend!, trajectory!, trajectory_hist!, plot

include(joinpath(@__DIR__, "..", "plotting", "plotting-itp.jl"))
export check_land, check_itp

include("io.jl")
export rtr2mat, rtr2nc

include("show.jl")
export length, show, iterate # various Base extensions

# initialize interpolants
function __init__()
    include(joinpath(@__DIR__, "..", "interpolants", "data", "data.jl"))
    itps_default_construct()
    itps_default_assign()
end

end # module

