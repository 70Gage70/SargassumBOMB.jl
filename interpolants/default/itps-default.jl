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

function _dtk_path(name::String)
    try
        open(dataset(name), DataToolkit.FilePath).path
    catch
        nothing
    end
end

"""
    itps_default_construct(; download_data = false, verbose = true)

Construct all interpolants using the default data. This overwrites any default interpolants already constructed.

Interpolants constructed: water, wind, stokes, waves, nutrients, temperature, land.

### Optional Arguments

- `download_data`: If `true`, the data required to construct the interpolants will be downloaded using `DataToolkit.jl`. \
This is roughly 1.2 GB of .nc files. Default `false`.
- `verbose`: If `true`, print itp construction stats. Default `true`.
"""
function itps_default_construct(; download_data::Bool = false, verbose::Bool = true)
    path2datatoml = joinpath(@__DIR__, "Data.toml") |> abspath
    data_col = loadcollection!(path2datatoml, @__MODULE__)


    if download_data
        DataToolkitCommon.Store.fetch!(data_col) # v0.9
        # DataToolkitStore.fetch!(data_col) # v0.10
    end

    verbose && @info "Constructing default interpolants."

    missings = String[]

    ########################### WATER
    verbose && @info "Constructing WATER interpolant."
    infile = _dtk_path("WATER")
    if infile !== nothing
        outfile = joinpath(@__DIR__, "itps", "WATER_ITP.jld2")
        rm(outfile, force = true)
        _construct_water_default(infile, outfile)
    else
        push!(missings, "WATER")
    end

    ########################### WIND
    verbose && @info "Constructing WIND interpolant."
    infile = _dtk_path("WIND")
    if infile !== nothing
        outfile = joinpath(@__DIR__, "itps", "WIND_ITP.jld2")
        rm(outfile, force = true)
        _construct_wind_default(infile, outfile)
    else
        push!(missings, "WIND")
    end

    ########################### STOKES
    verbose && @info "Constructing STOKES interpolant."
    infile = _dtk_path("WAVES")
    if infile !== nothing
        outfile = joinpath(@__DIR__, "itps", "STOKES_ITP.jld2")
        rm(outfile, force = true)
        _construct_stokes_default(infile, outfile)
    else
        push!(missings, "STOKES")
    end

    ########################### WAVES
    verbose && @info "Constructing WAVES interpolant."
    infile = _dtk_path("WAVES")
    if infile !== nothing
        outfile = joinpath(@__DIR__, "itps", "WAVES_ITP.jld2")
        rm(outfile, force = true)
        _construct_waves_default(infile, outfile)
    else
        push!(missings, "WAVES")
    end

    ########################### NUTRIENTS
    verbose && @info "Constructing NUTRIENTS interpolant."
    infile = _dtk_path("NUTRIENTS")
    if infile !== nothing
        outfile = joinpath(@__DIR__, "itps", "NUTRIENTS_ITP.jld2")
        rm(outfile, force = true)
        _construct_nutrients_default(infile, outfile)
    else
        push!(missings, "NUTRIENTS")
    end

    ########################### TEMPERATURE
    verbose && @info "Constructing TEMPERATURE interpolant."
    infile = _dtk_path("TEMPERATURE")
    if infile !== nothing
        outfile = joinpath(@__DIR__, "itps", "TEMPERATURE_ITP.jld2")
        rm(outfile, force = true)
        _construct_temperature_default(infile, outfile)
    else
        push!(missings, "TEMPERATURE")
    end

    ########################### LAND
    verbose && @info "Constructing LAND interpolant."
    try 
        outfile = joinpath(@__DIR__, "itps", "LAND_ITP.jld2")
        rm(outfile, force = true)
        _construct_land_default(outfile)
    catch
        @warn "Could not construct land interpolant."
    end

    ########################### END MATTER

    if length(missings) == 6
        verbose && @warn "Could not construct any interpolants; data missing. Try running `itps_default_construct(download_data = true)`. This downloads roughly 1.2 GB of data."
    elseif 0 < length(missings) < 6
        verbose && @warn "Could not construct interpolants $(missings); data missing. Try running `itps_default_construct(download_data = true)`. This downloads roughly 1.2 GB of data."
    else
        itps_load(ITPS_DEFAULT_DIR)
        verbose && @info "Default interpolants constructed."
    end

    return nothing
end