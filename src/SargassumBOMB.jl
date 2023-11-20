module SargassumBOMB

# core functionality
using OrdinaryDiffEq, SargassumFromAFAI, NearestNeighbors, Interpolations
using Unitful, Surrogates, Dates
using LinearAlgebra: norm, ⋅

# i/o
using MAT, NetCDF, JLD2

# probability/statistics
using StatsBase, Distributions 
using Random: seed!

# plotting/geography
using Makie, CairoMakie, GeoMakie
using GeoMakie, GeoMakie.GeoJSON, GeoDatasets, GeoInterface

# printing
using Crayons.Box
using Latexify

include("coordinates.jl")
export EquirectangularReference, EQR_DEFAULT, sph2xy, xy2sph

include(joinpath(@__DIR__, "..", "interpolants", "itp-core.jl"))
export GriddedField, InterpolatedField, interpolate, add_derivatives, reduce_vector_to_range, rata2datetime_minute

include(joinpath(@__DIR__, "..", "interpolants", "ITPConstruct.jl"))
export construct_all_itp
export WATER_ITP, WIND_ITP, WAVES_ITP, STOKES_ITP, LAND_ITP, TEMP_ITP, NO3_ITP

include(joinpath(@__DIR__, "biology.jl"))
export AbstractGrowthDeathModel, ImmortalModel, BrooksModelParameters, brooks_dSdt_clump, brooks_dSdt_raft, BrooksModel

include("geography.jl")
export AbstractLand, NoLand, Land

include("raft-parameters.jl")
export ClumpParameters
export SpringParameters, spring_force 
export InitialConditions
export AbstractConnections, ConnectionsNone, ConnectionsFull, ConnectionsRadius, ConnectionsNearest, form_connections!
export RaftParameters

include("physics.jl")
export Raft!, WaterWind!

include("control.jl")
export n_clumps, clump_i, com, cb_update, cb_land, cb_growth_death, cb_connections, kill!, grow!

include("trajectories.jl")
export Trajectory, time_slice, RaftTrajectory, uniformize, bins

include(joinpath(@__DIR__, "..", "plotting", "plotting-core.jl"))
export default_fig, geo_axis, land!, data_legend!, trajectory!, vector_field_t!, scalar_field_t!, trajectory_hist!

include(joinpath(@__DIR__, "..", "plotting", "plotting-itp.jl"))
export check_land, check_windwater

include("main.jl")
export simulate

export length, show, iterate # various Base extensions

# initialize interpolants
function __init__()
    construct_all_itp()
end

import PrecompileTools

PrecompileTools.@compile_workload begin
    construct_all_itp()
end

PrecompileTools.@compile_workload begin
    ics = InitialConditions(range(-55.0, -50.0, length = 5), range(5.0, 10.0, length = 5), ref = EQR_DEFAULT)
    clumps = ClumpParameters(EQR_DEFAULT)
    springs = SpringParameters(k -> 0.1, 100.0)
    connections = ConnectionsNearest(10)
    gd_model = BrooksModel()
    land = Land()
    tspan = (0.0, 5.0)

    rp = RaftParameters(;
    ics = ics,
    clumps = clumps,
    springs = springs,
    connections = connections,
    t0 = first(tspan),
    gd_model = gd_model,
    land = land)

    sol = simulate(rp, tspan, showprogress = false)
end

end # module

