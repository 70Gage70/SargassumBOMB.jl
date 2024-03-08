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
    const STOKES_ITP

The interpolant for Stokes drift velocity. This is a `Ref`, access or modify the actual interpolant with `STOKES_ITP.x`.
"""
const STOKES_ITP = Ref{InterpolatedField}()

"""
    const WAVES_ITP

The interpolant for wave height. This is a `Ref`, access or modify the actual interpolant with `WAVES_ITP.x`.
"""
const WAVES_ITP = Ref{InterpolatedField}()

"""
    const NUTRIENTS_ITP

The interpolant for ocean nitrogen content. This is a `Ref`, access or modify the actual interpolant with `NO3_ITP_ITP.x`.
"""
const NUTRIENTS_ITP = Ref{InterpolatedField}()

"""
    const TEMPERATURE_ITP

The interpolant for ocean temperature. This is a `Ref`, access or modify the actual interpolant with `TEMPERATURE_ITP.x`.
"""
const TEMPERATURE_ITP = Ref{InterpolatedField}()

"""
    const LAND_ITP

The interpolant for landmass location. This is a `Ref`, access or modify the actual interpolant with `LAND_ITP.x`.
"""
const LAND_ITP = Ref{InterpolatedField}()


"""
    itps_default_construct(; force = false)

Construct all interpolants using the default data.

Create for the file `~/interpolants/itps/itps-constructed.jl` to flag that the construction has been completed previously. If 
it has, no action is performed. Bypass this by passing `force == true`, in which case all interpolants are recreated from scratch.

Interpolants constructed: water, wind, stokes, waves, nutrients, temperature, land.
"""
function itps_default_construct(; force::Bool = false)
    itp_constructed_flag = joinpath(@__DIR__, "itps", "itps-constructed.jl")
    if !force && isfile(itp_constructed_flag)
        return nothing
    end

    @info "Constructing default interpolants."

    ###########################
    ########################### WATER
    ###########################

    infile = joinpath(@__DIR__, "data", "water.nc")
    outfile = joinpath(@__DIR__, "itps", "WATER_ITP.jld2")

    println("INFILE")
    println(infile)

    rm(outfile, force = true)

    gf = GriddedField(3)

    add_spatial_dimension!(gf, infile, "lon", :lon, u"°", "degrees")
    add_spatial_dimension!(gf, infile, "lat", :lat, u"°", "degrees")
    add_temporal_dimension!(gf, infile, "time", :t, DateTime(0000, 12, 31), Day) # Rata Die days

    add_field!(gf, infile, "u", :u, u"km/d", "speed")
    add_field!(gf, infile, "v", :v, u"km/d", "speed")

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf)
    add_derivatives!(itp)

    jldsave(outfile, WATER_ITP = itp)

    ###########################
    ########################### WIND
    ###########################

    infile = joinpath(@__DIR__, "data", "wind.nc")
    outfile = joinpath(@__DIR__, "itps", "WIND_ITP.jld2")
    rm(outfile, force = true)

    gf = GriddedField(3)

    add_spatial_dimension!(gf, infile, "longitude", :lon, u"°", "degrees")
    add_spatial_dimension!(gf, infile, "latitude", :lat, u"°", "degrees")
    add_temporal_dimension!(gf, infile, "time", :t, DateTime(1900, 1, 1), Hour)

    add_field!(gf, infile, "u10", :u, u"m/s", "speed")
    add_field!(gf, infile, "v10", :v, u"m/s", "speed")

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf)
    add_derivatives!(itp)

    jldsave(outfile, WIND_ITP = itp)

    ###########################
    ########################### STOKES
    ###########################

    infile = joinpath(@__DIR__, "data", "waves.nc")
    outfile = joinpath(@__DIR__, "itps", "STOKES_ITP.jld2")
    rm(outfile, force = true)

    gf = GriddedField(3)

    add_spatial_dimension!(gf, infile, "longitude", :lon, u"°", "degrees")
    add_spatial_dimension!(gf, infile, "latitude", :lat, u"°", "degrees")
    add_temporal_dimension!(gf, infile, "time", :t, DateTime(1900, 1, 1), Hour)

    add_field!(gf, infile, "ust", :u, u"m/s", "speed")
    add_field!(gf, infile, "vst", :v, u"m/s", "speed")

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf)
    add_derivatives!(itp)

    jldsave(outfile, STOKES_ITP = itp)

    ###########################
    ########################### WAVES
    ###########################

    infile = joinpath(@__DIR__, "data", "waves.nc")
    outfile = joinpath(@__DIR__, "itps", "WAVES_ITP.jld2")
    rm(outfile, force = true)

    gf = GriddedField(3)

    add_spatial_dimension!(gf, infile, "longitude", :lon, u"°", "degrees")
    add_spatial_dimension!(gf, infile, "latitude", :lat, u"°", "degrees")
    add_temporal_dimension!(gf, infile, "time", :t, DateTime(1900, 1, 1), Hour)

    add_field!(gf, infile, "swh", :swh, u"m", "wave_height")

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf)
    # add_derivatives!(itp)

    jldsave(outfile, WAVES_ITP = itp)

    ###########################
    ########################### NUTRIENTS
    ###########################

    infile = joinpath(@__DIR__, "data", "nutrients.nc")
    outfile = joinpath(@__DIR__, "itps", "NUTRIENTS_ITP.jld2")
    rm(outfile, force = true)

    gf = GriddedField(3)

    add_spatial_dimension!(gf, infile, "longitude", :lon, u"°", "degrees")
    add_spatial_dimension!(gf, infile, "latitude", :lat, u"°", "degrees")
    add_temporal_dimension!(gf, infile, "time", :t, DateTime(1950, 1, 1), Hour)

    add_field!(gf, infile, "no3", :no3, u"mmol/m^3", "concentration", take_axes = [:,:,1,:])

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf)
    # add_derivatives!(itp)

    jldsave(outfile, NUTRIENTS_ITP = itp)

    ###########################
    ########################### TEMPERATURE
    ###########################

    infile = joinpath(@__DIR__, "data", "currents-temperature.nc")
    outfile = joinpath(@__DIR__, "itps", "TEMPERATURE_ITP.jld2")
    rm(outfile, force = true)

    gf = GriddedField(3)

    add_spatial_dimension!(gf, infile, "longitude", :lon, u"°", "degrees")
    add_spatial_dimension!(gf, infile, "latitude", :lat, u"°", "degrees")
    add_temporal_dimension!(gf, infile, "time", :t, DateTime(1950, 1, 1), Day)

    add_field!(gf, infile, "thetao_glor", :temp, u"°C", "temperature", take_axes = [:,:,1,:])

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf)
    # add_derivatives!(itp)

    jldsave(outfile, TEMPERATURE_ITP = itp)

    ###########################
    ########################### LAND
    ###########################

    outfile = joinpath(@__DIR__, "itps", "LAND_ITP.jld2")
    rm(outfile, force = true)

    gf = GriddedField(2)

    # 0 is ocean, 1 is land and 2 is lake
    lon, lat, land = GeoDatasets.landseamask(resolution = 'c', grid = 5)
    land[land .== 2] .= 1 # lake is not ocean, so it's land

    push!(gf.dims_names, (:lon, u"°"))
    gf.dims[:lon] = lon
    push!(gf.dims_names, (:lat, u"°"))
    gf.dims[:lat] = lat

    push!(gf.fields_names, (:land, NoUnits))
    gf.fields[:land] = land

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf, interpolant_type = "nearest")

    jldsave(outfile, LAND_ITP = itp)

    ###########################
    ########################### END MATTER
    ###########################

    touch(itp_constructed_flag)

    @info "Default interpolants constructed."

    return nothing
end

"""
    itps_default_assign(; force = false)

Attempt to assign each interpolant to its default value, constructed from [`itp_default_construct`](@ref).
"""
function itps_default_assign()
    try
        WATER_ITP.x = load(joinpath(@__DIR__, "itps", "WATER_ITP.jld2"), "WATER_ITP")
    catch
        @warn "Could not assign `WATER_ITP`."
    end

    try
        WIND_ITP.x = load(joinpath(@__DIR__, "itps", "WIND_ITP.jld2"), "WIND_ITP")
    catch
        @warn "Could not assign `WIND_ITP`."
    end

    try
        STOKES_ITP.x = load(joinpath(@__DIR__, "itps", "STOKES_ITP.jld2"), "STOKES_ITP")
    catch
        @warn "Could not assign `STOKES_ITP`."
    end

    try
        WAVES_ITP.x = load(joinpath(@__DIR__, "itps", "WAVES_ITP.jld2"), "WAVES_ITP")
    catch
        @warn "Could not assign `WAVES_ITP`."
    end

    try
        NUTRIENTS_ITP.x = load(joinpath(@__DIR__, "itps", "NUTRIENTS_ITP.jld2"), "NUTRIENTS_ITP")
    catch
        @warn "Could not assign `NUTRIENTS_ITP`."
    end

    try
        TEMPERATURE_ITP.x = load(joinpath(@__DIR__, "itps", "TEMPERATURE_ITP.jld2"), "TEMPERATURE_ITP")
    catch
        @warn "Could not assign `TEMPERATURE_ITP`."
    end

    try
        LAND_ITP.x = load(joinpath(@__DIR__, "itps", "LAND_ITP.jld2"), "LAND_ITP")
    catch
        @warn "Could not assign `LAND_ITP`."
    end

    return nothing
end