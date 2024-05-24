@RemoteFile(
    SBOMB_WATER, 
    "https://www.dropbox.com/scl/fi/7lvc5j0dqcbztcwd3rnnp/water.nc?rlkey=yhpmf5evtzrppr14v691l97j3&dl=1",
    file = "water.nc", 
    updates=:never,
    dir = joinpath(@__DIR__, "data"))

@RemoteFile(
    SBOMB_WIND, 
    "https://www.dropbox.com/scl/fi/k3vrbuxegp6mj6y5oqpbo/wind.nc?rlkey=jgi2zvrwagc2e1tcwt1nvzz54&dl=0",
    file = "wind.nc", 
    updates=:never,
    dir = joinpath(@__DIR__, "data"))

@RemoteFile(
    SBOMB_WAVES, 
    "https://www.dropbox.com/scl/fi/j5b6akflmgvhcfy124l4l/waves.nc?rlkey=me7c7tne70kd5xrt3oy7o15tm&dl=0",
    file = "waves.nc", 
    updates=:never,
    dir = joinpath(@__DIR__, "data"))

@RemoteFile(
    SBOMB_TEMP, 
    "https://www.dropbox.com/scl/fi/oyrk1klf8rugk3mty0erx/currents-temperature.nc?rlkey=3t6znxletyret5jo0v8447dpf&dl=1",
    file = "currents-temperature.nc", 
    updates=:never,
    dir = joinpath(@__DIR__, "data"))

@RemoteFile(
    SBOMB_NUTRIENTS, 
    "https://www.dropbox.com/scl/fi/pj3ovx8owlxqxaj43x2py/nutrients.nc?rlkey=jsulyltckanjqdmi8ygxrddwo&dl=1",
    file = "nutrients.nc", 
    updates=:never,
    dir = joinpath(@__DIR__, "data"))


function _construct_water_default(infile::String, outfile::String)
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

    return nothing
end

function _construct_wind_default(infile::String, outfile::String)
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

    return nothing
end

function _construct_stokes_default(infile::String, outfile::String)
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

    return nothing
end

function _construct_waves_default(infile::String, outfile::String)
    gf = GriddedField(3)

    add_spatial_dimension!(gf, infile, "longitude", :lon, u"°", "degrees")
    add_spatial_dimension!(gf, infile, "latitude", :lat, u"°", "degrees")
    add_temporal_dimension!(gf, infile, "time", :t, DateTime(1900, 1, 1), Hour)

    add_field!(gf, infile, "swh", :swh, u"m", "wave_height")

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf, interpolant_type = "nearest")

    jldsave(outfile, WAVES_ITP = itp)

    return nothing
end

function _construct_nutrients_default(infile::String, outfile::String)
    gf = GriddedField(3)

    add_spatial_dimension!(gf, infile, "longitude", :lon, u"°", "degrees")
    add_spatial_dimension!(gf, infile, "latitude", :lat, u"°", "degrees")
    add_temporal_dimension!(gf, infile, "time", :t, DateTime(1950, 1, 1), Hour)

    add_field!(gf, infile, "no3", :no3, u"mmol/m^3", "concentration", take_axes = [:,:,1,:])

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf, interpolant_type = "nearest")

    jldsave(outfile, NUTRIENTS_ITP = itp)

    return nothing
end

function _construct_temperature_default(infile::String, outfile::String)
    gf = GriddedField(3)

    add_spatial_dimension!(gf, infile, "longitude", :lon, u"°", "degrees")
    add_spatial_dimension!(gf, infile, "latitude", :lat, u"°", "degrees")
    add_temporal_dimension!(gf, infile, "time", :t, DateTime(1950, 1, 1), Day)

    add_field!(gf, infile, "thetao_glor", :temp, u"°C", "temperature", take_axes = [:,:,1,:])

    ranges_increasing!(gf)
    sph2xy!(gf)

    itp = InterpolatedField(gf, interpolant_type = "nearest")

    jldsave(outfile, TEMPERATURE_ITP = itp)

    return nothing
end

function _construct_land_default(outfile::String)
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

    return nothing
end


"""
    itps_default_construct(; download_data = false)

Construct all interpolants using the default data. This overwrites any default interpolants already constructed.

Interpolants constructed: water, wind, stokes, waves, nutrients, temperature, land.

### Optional Arguments

- `download_data`: If `true`, the data required to construct the interpolants will be downloaded using `RemoteFiles.jl`. \
This is roughly 1 GB of .nc files. Any existing data will be overwritten. Default `false`.
"""
function itps_default_construct(; download_data::Bool = false)
    if download_data
        download(SBOMB_WATER)
        download(SBOMB_WIND)
        download(SBOMB_WAVES)
        download(SBOMB_TEMP)
        download(SBOMB_NUTRIENTS)
    end

    @info "Constructing default interpolants."

    missings = String[]

    ########################### WATER
    infile = joinpath(@__DIR__, "data", "water.nc")
    if isfile(infile)
        outfile = joinpath(@__DIR__, "itps", "WATER_ITP.jld2")
        rm(outfile, force = true)
        _construct_water_default(infile, outfile)
    else
        push!(missings, "WATER")
    end

    ########################### WIND
    infile = joinpath(@__DIR__, "data", "wind.nc")
    if isfile(infile)
        outfile = joinpath(@__DIR__, "itps", "WIND_ITP.jld2")
        rm(outfile, force = true)
        _construct_wind_default(infile, outfile)
    else
        push!(missings, "WIND")
    end

    ########################### STOKES
    infile = joinpath(@__DIR__, "data", "waves.nc")
    if isfile(infile)
        outfile = joinpath(@__DIR__, "itps", "STOKES_ITP.jld2")
        rm(outfile, force = true)
        _construct_stokes_default(infile, outfile)
    else
        push!(missings, "STOKES")
    end

    ########################### WAVES
    infile = joinpath(@__DIR__, "data", "waves.nc")
    if isfile(infile)
        outfile = joinpath(@__DIR__, "itps", "WAVES_ITP.jld2")
        rm(outfile, force = true)
        _construct_waves_default(infile, outfile)
    else
        push!(missings, "WAVES")
    end

    ########################### NUTRIENTS
    infile = joinpath(@__DIR__, "data", "nutrients.nc")
    if isfile(infile)
        outfile = joinpath(@__DIR__, "itps", "NUTRIENTS_ITP.jld2")
        rm(outfile, force = true)
        _construct_nutrients_default(infile, outfile)
    else
        push!(missings, "NUTRIENTS")
    end

    ########################### TEMPERATURE
    infile = joinpath(@__DIR__, "data", "currents-temperature.nc")
    if isfile(infile)
        outfile = joinpath(@__DIR__, "itps", "TEMPERATURE_ITP.jld2")
        rm(outfile, force = true)
        _construct_temperature_default(infile, outfile)
    else
        push!(missings, "TEMPERATURE")
    end

    ########################### LAND
    try 
        outfile = joinpath(@__DIR__, "itps", "LAND_ITP.jld2")
        rm(outfile, force = true)
        _construct_land_default(outfile)
    catch
        @warn "Could not construct land interpolant."
    end

    ########################### END MATTER

    if length(missings) == 6
        @warn "Could not construct any interpolants; data missing. Try running `itps_default_construct(download_data = true)`"
    elseif 0 < length(missings) < 6
        @warn "Could not construct interpolants $(missings); data missing. Try running `itps_default_construct(download_data = true)`"
    else
        itps_load(ITPS_DEFAULT_DIR)
        @info "Default interpolants constructed."
    end

    return nothing
end