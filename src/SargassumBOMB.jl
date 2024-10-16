module SargassumBOMB

# core functionality
using OrdinaryDiffEqTsit5, SargassumFromAFAI, NearestNeighbors, Interpolations
using SparseArrays
using Unitful, Dates
using LinearAlgebra: norm, ⋅ # dot product
using ArgCheck

# i/o
using MAT, NetCDF, JLD2

# downloading interpolants
using Scratch, Downloads
const _ITPS_RAW_SCRATCH = Ref{String}()
const _ITPS_SCRATCH = Ref{String}()

# probability/statistics
using StatsBase, Distributions 
using Random: seed!

# plotting/geography
using Makie, CairoMakie
using SargassumColors
using GeoDatasets # landseamask

# printing
using Crayons.Box
using Latexify
using ProgressBars

# optimization 
using Metaheuristics
using QuasiMonteCarlo

###############################################################

include("utils.jl")
export n_clumps, clump_i, com, vec2range

include("coordinates.jl")
export UNITS, EARTH_RADIUS, EquirectangularReference, EQR
export sph2xy, xy2sph, γ_sphere, τ_sphere

include("time.jl")
export T_REF, datetime2time, time2datetime, ymw2time, time2ymw, ymwspan2weekspan, ymwplusweek

include(joinpath(@__DIR__, "..", "interpolants", "itps-core.jl"))
export GriddedField, InterpolatedField
export add_spatial_dimension!, add_temporal_dimension!, add_field!, ranges_increasing!, sph2xy!, add_derivatives!

include(joinpath(@__DIR__, "..", "interpolants", "itps-definitions.jl"))
export ITPS_DEFAULT_DIR
export WATER_ITP, WIND_ITP, WAVES_ITP, STOKES_ITP, LAND_ITP, TEMPERATURE_ITP, NUTRIENTS_ITP
export itps_load

include(joinpath(@__DIR__, "..", "interpolants", "default", "itps-default.jl"))
export itps_default_construct

include(joinpath(@__DIR__, "..", "interpolants", "itps-interface.jl"))
export update_interpolant!, limits, dim, dims, field, fields

include("land.jl")
export AbstractLand, NoLand, Land

include("ics.jl")
export InitialConditions

include("springs.jl")
export AbstractSpring, HookeSpring, BOMBSpring, ΔL, spring_force
export AbstractConnections, ConnectionsNone, ConnectionsFull, ConnectionsRadius, ConnectionsNearest, form_connections

include("growth-death.jl")
export AbstractGrowthDeathModel, ImmortalModel, BrooksModelParameters, BrooksModel

include("rafts-clumps.jl")
export ClumpParameters, RaftParameters, dxdy_MR

include("physics.jl")
export FastRaft!, Raft!, Leeway!

include("control.jl")
export kill!, grow!

include("trajectories.jl")
export Trajectory, time_slice, RaftTrajectory, bins

include("main.jl")
export simulate

include("optimization.jl")
export TimeSeries, vec
export LossFunction
export optimize!

include("plotting-core.jl")
export trajectory!, trajectory_hist!, trajectory

include("plotting-itp.jl")
export check_land, check_itp

include("io.jl")
export rtr2mat, rtr2nc

include("show.jl")
export length, show, iterate # various Base extensions

include(joinpath(@__DIR__, "..", "examples", "examples.jl"))
export Examples

# initialize interpolants
function __init__()
    _ITPS_RAW_SCRATCH.x = @get_scratch!("_ITPS_RAW_SCRATCH")
    _ITPS_SCRATCH.x = @get_scratch!("_ITPS_SCRATCH")

    try
        itps_load()
    catch
        try 
            itps_default_construct()
        catch
            nothing
        end
    end

    return nothing
end

end # module

