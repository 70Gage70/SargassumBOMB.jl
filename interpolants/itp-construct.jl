# module ITPConstruct

# using LinearAlgebra: norm, ⋅
# using MAT
# using Interpolations
# using Dates
# using JLD2
# using NetCDF
# using Unitful

include(joinpath(@__DIR__, "..", "data", "preprocessed", "preprocess-all.jl"))

include(joinpath(@__DIR__, "biology", "itp-construct.jl"))
include(joinpath(@__DIR__, "land", "itp-construct.jl"))
include(joinpath(@__DIR__, "ocean-atmos", "itp-construct.jl"))
include(joinpath(@__DIR__, "waves", "itp-construct.jl"))

"""
    const WATER_ITP

The interpolant for ocean currents. This is a `Ref`, access or modify the actual interpolant with `WATER_ITP.x`.
"""
const WATER_ITP = Ref{InterpolatedField}()

"""
    const WIND_ITP

The interpolant for wind speed. This is a `Ref`, access or modify the actual interpolant with `WIND_ITP.x`.
"""
const WIND_ITP = Ref{InterpolatedField}()

"""
    const WAVES_ITP

The interpolant for wave height. This is a `Ref`, access or modify the actual interpolant with `WAVES_ITP.x`.
"""
const WAVES_ITP = Ref{InterpolatedField}()

"""
    const STOKES_ITP

The interpolant for Stokes drift velocity. This is a `Ref`, access or modify the actual interpolant with `STOKES_ITP.x`.
"""
const STOKES_ITP = Ref{InterpolatedField}()

"""
    const LAND_ITP

The interpolant for landmass location. This is a `Ref`, access or modify the actual interpolant with `LAND_ITP.x`.
"""
const LAND_ITP = Ref{InterpolatedField}()

"""
    const TEMP_ITP

The interpolant for ocean temperature. This is a `Ref`, access or modify the actual interpolant with `TEMP_ITP.x`.
"""
const TEMP_ITP = Ref{InterpolatedField}()

"""
    const NO3_ITP

The interpolant for ocean nitrogen content. This is a `Ref`, access or modify the actual interpolant with `NO3_ITP_ITP.x`.
"""
const NO3_ITP = Ref{InterpolatedField}()

"""
    construct_all_itp(;force_preprocess = false, force_itp = false)

Preprocess the raw data using [`preprocess_all`](@ref) and construct all interpolants.

### Optional Arguments

If `force_preprocess == true`, the data will be preprocessed even if preprocessed data already exists. Default `false`.

If `force == true`, the interpolants will be computed even if they already exist. Default `false`.
"""
function construct_all_itp(;force_preprocess::Bool = false, force_itp::Bool = false)
    preprocess_path = joinpath(@__DIR__, "..", "data", "preprocessed", "is-preprocessed.jl")
    if !isfile(preprocess_path) || force_preprocess
        preprocess_all()
        touch(preprocess_path)
    end

    water_path = joinpath(@__DIR__, "..", "interpolants", "ocean-atmos", "WATER_ITP.jld2")
    if !isfile(water_path) || force_itp
        construct_water_itp()
    end
    WATER_ITP.x = load(water_path, "WATER_ITP")

    wind_path = joinpath(@__DIR__, "..", "interpolants", "ocean-atmos", "WIND_ITP.jld2")
    if !isfile(wind_path) || force_itp
        construct_wind_itp()
    end
    WIND_ITP.x = load(wind_path, "WIND_ITP")

    waves_path = joinpath(@__DIR__, "..", "interpolants", "waves", "WAVES_ITP.jld2")
    if !isfile(waves_path) || force_itp
        construct_waves_itp()
    end
    WAVES_ITP.x = load(waves_path, "WAVES_ITP")

    stokes_path = joinpath(@__DIR__, "..", "interpolants", "waves", "STOKES_ITP.jld2")
    if !isfile(stokes_path) || force_itp
        construct_stokes_itp()
    end
    STOKES_ITP.x = load(stokes_path, "STOKES_ITP")

    land_path = joinpath(@__DIR__, "..", "interpolants", "land", "LAND_ITP.jld2")
    if !isfile(land_path) || force_itp
        construct_land_itp()
    end
    LAND_ITP.x = load(land_path, "LAND_ITP")

    temp_path = joinpath(@__DIR__, "..", "interpolants", "biology", "TEMP_ITP.jld2")
    if !isfile(temp_path) || force_itp
        construct_temp_itp() 
    end
    TEMP_ITP.x = load(temp_path, "TEMP_ITP")

    no3_path = joinpath(@__DIR__, "..", "interpolants", "biology", "NO3_ITP.jld2")
    if !isfile(no3_path) || force_itp
        construct_no3_itp()
    end
    NO3_ITP.x = load(no3_path, "NO3_ITP")

    return nothing
end



# export construct_all_itp
# export 

# end # module
