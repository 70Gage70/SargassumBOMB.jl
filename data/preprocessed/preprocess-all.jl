const UNITS_OUT_HEIGHT = u"m"
const UNITS_OUT_SPEED = u"km/d"
const UNITS_OUT_TEMP = u"°C"
const UNITS_OUT_NUTR = u"mmol/m^3"

const T_START_ITP = DateTime(2018, 1, 1)
const T_END_ITP = DateTime(2018, 12, 31)

"""
    preprocess_all()

Compute all the preprocessed `.mat` files from raw data. Any preexisting `.mat` files will be replaced.

Generally, [`construct_all_itp`](@ref) should be used instead of this function.
"""
function preprocess_all()
    raw_path = joinpath(@__DIR__, "..", "raw")

    wind_raw = joinpath(raw_path, "wind.nc")
    water_temp_raw = joinpath(raw_path, "currents-temperature.nc")
    nutr_raw = joinpath(raw_path, "nutrients.nc")
    water_raw = joinpath(raw_path, "water.nc")
    waves_raw = joinpath(raw_path, "waves.nc")

    NetCDF_units = Dict(
        "m" => u"m",
        "m s**-1" => u"m/s",
        "m s-1" => u"m/s", 
        "degrees_C" => u"°C",
        "mmol m-3" => u"mmol/m^3",
        "km/d" => u"km/d"
    )

    # use ncinfo(file) to see attributes/metadata
    # use ncgetatt(file, variable, attribute) to read an attribute of a variable
    # use ncread(file, variable) to read the value of a variable

    ###########################
    ########################### WIND
    ###########################

    # variables: longitude, latitude, time, u10(longitude, latitude, time), v10(longitude, latitude, time)
    # longitude is in degrees E/W [-180, 180]
    # latitude is in degrees N/S [-90, 90]
    # time is in hours since 1900-01-01 00:00:00
    # u10 and v10 are in m/s

    lon_wind = ncread(wind_raw, "longitude") .|> Float64
    lat_wind = ncread(wind_raw, "latitude") .|> Float64

    tref = DateTime(1900, 1, 1, 0, 0, 0)
    time_wind = ncread(wind_raw, "time") .|> (x -> tref + Hour(x)) .|> datetime2rata
    T_START_ITP_idx = searchsortedfirst(rata2datetime.(time_wind), T_START_ITP)
    T_END_ITP_idx = searchsortedfirst(rata2datetime.(time_wind), T_END_ITP)
    time_wind = time_wind[T_START_ITP_idx:T_END_ITP_idx]

    # u

    u_wind = ncread(wind_raw, "u10")[:, :, T_START_ITP_idx:T_END_ITP_idx]
    scale_factor = ncgetatt(wind_raw, "u10", "scale_factor")
    units_factor = ncgetatt(wind_raw, "u10", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
    fill_value = ncgetatt(wind_raw, "u10", "_FillValue")

    u_wind[u_wind .== fill_value] .= 0.0
    u_wind = u_wind .|> x -> scale_factor * units_factor * x 

    # v

    v_wind = ncread(wind_raw, "v10")[:, :, T_START_ITP_idx:T_END_ITP_idx]
    scale_factor = ncgetatt(wind_raw, "v10", "scale_factor")
    units_factor = ncgetatt(wind_raw, "v10", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
    fill_value = ncgetatt(wind_raw, "v10", "_FillValue")

    v_wind[v_wind .== fill_value] .= 0.0
    v_wind = v_wind .|> x -> scale_factor * units_factor * x 

    ###########################
    ########################### TEMPERATURE
    ###########################

    # variables: longitude, latitude, depth, time, thetao_glor(longitude, latitude, depth, time)
    # longitude is in degrees E/W [-180, 180]
    # latitude is in degrees N/S [-90 - 90]
    # depth is in m
    # time is in days since 1950-01-01 00:00:00
    # thetao_glor is in degrees_C

    lon_temp = ncread(water_temp_raw, "longitude") .|> Float64 
    lat_temp = ncread(water_temp_raw, "latitude") .|> Float64

    tref = DateTime(1950, 1, 1, 0, 0, 0)
    time_temp = ncread(water_temp_raw, "time") .|> (x -> tref + Day(x)) .|> datetime2rata
    T_START_ITP_idx = searchsortedfirst(rata2datetime.(time_temp), T_START_ITP)
    T_END_ITP_idx = searchsortedfirst(rata2datetime.(time_temp), T_END_ITP)
    time_temp = time_temp[T_START_ITP_idx:T_END_ITP_idx]

    # temp

    temp = ncread(water_temp_raw, "thetao_glor")[:, :, 1, T_START_ITP_idx:T_END_ITP_idx]
    scale_factor = ncgetatt(water_temp_raw, "thetao_glor", "scale_factor")
    units_factor = ncgetatt(water_temp_raw, "thetao_glor", "units") |> x -> uconvert(UNITS_OUT_TEMP, 1.0 * NetCDF_units[x]).val
    fill_value = ncgetatt(water_temp_raw, "thetao_glor", "_FillValue")

    temp[temp .== fill_value] .= 0.0
    temp = temp .|> x -> scale_factor * units_factor * x 

    ###########################
    ########################### NUTRIENTS
    ###########################

    # variables: longitude, latitude, depth, time, no3(longitude, latitude, depth, time)
    # longitude is in degrees E/W [-180, 180]
    # latitude is in degrees N/S [-90 - 90]
    # depth is in m
    # time is in hours since 1950-01-01 00:00:00
    # no3 is in mmol/m^3

    lon_nutr = ncread(nutr_raw, "longitude") .|> Float64 
    lat_nutr = ncread(nutr_raw, "latitude") .|> Float64

    tref = DateTime(1950, 1, 1, 0, 0, 0)
    time_nutr = ncread(nutr_raw, "time") .|> (x -> tref + Hour(x)) .|> datetime2rata
    T_START_ITP_idx = searchsortedfirst(rata2datetime.(time_nutr), T_START_ITP)
    T_END_ITP_idx = searchsortedfirst(rata2datetime.(time_nutr), T_END_ITP)
    time_nutr = time_nutr[T_START_ITP_idx:T_END_ITP_idx]

    # no3

    no3 = ncread(nutr_raw, "no3")[:, :, 1, T_START_ITP_idx:T_END_ITP_idx]
    scale_factor = 1.0
    units_factor = ncgetatt(nutr_raw, "no3", "units") |> x -> uconvert(UNITS_OUT_NUTR, 1.0 * NetCDF_units[x]).val
    fill_value = ncgetatt(nutr_raw, "no3", "_FillValue")

    no3[no3 .== fill_value] .= 0.0
    no3 = no3 .|> x -> scale_factor * units_factor * x

    ###########################
    ########################### WATER
    ###########################

    # indep vars: "lon", "lat", "time" 
    # dep vars: "u", "v"
    # longitude is in degrees E/W [-180, 180]
    # latitude is in degrees N/S [-90 - 90]
    # time is in rata days
    # u/v are in km/d
    # contains no NaNs or scaling

    lon_water = ncread(water_raw, "lon") |> vec .|> Float64
    lon_idx = findall(x -> -101 <= x <= -39, lon_water)
    lon_water = lon_water[lon_idx]

    lat_water = ncread(water_raw, "lat") |> vec .|> Float64
    lat_idx = findall(x -> 0 <= x <= 40, lat_water)
    lat_water = lat_water[lat_idx]

    time_water = ncread(water_raw, "time") |> vec
    T_START_ITP_idx = searchsortedfirst(rata2datetime.(time_water), T_START_ITP)
    T_END_ITP_idx = searchsortedfirst(rata2datetime.(time_water), T_END_ITP)
    time_water = time_water[T_START_ITP_idx:T_END_ITP_idx]

    # u

    u_water = ncread(water_raw, "u")

    units_factor = ncgetatt(water_raw, "u", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
    u_water = u_water .|> x -> units_factor * x 

    # v

    v_water = ncread(water_raw, "v")

    units_factor = ncgetatt(water_raw, "v", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
    v_water = v_water .|> x -> units_factor * x 

    ###########################
    ########################### WAVES
    ###########################

    # variables: longitude, latitude, time, swh(longitude, latitude, time), ust(longitude, latitude, time), vst(longitude, latitude, time)
    # longitude is in degrees E/W [-180, 180]
    # latitude is in degrees N/S [-90, 90]
    # time is in hours since 1900-01-01 00:00:00
    # swh is in m
    # ust and vst are in m/s

    lon_waves = ncread(waves_raw, "longitude") .|> Float64
    lat_waves = ncread(waves_raw, "latitude") .|> Float64

    tref = DateTime(1900, 1, 1, 0, 0, 0)
    time_waves = ncread(waves_raw, "time") .|> (x -> tref + Hour(x)) .|> datetime2rata
    T_START_ITP_idx = searchsortedfirst(rata2datetime.(time_waves), T_START_ITP)
    T_END_ITP_idx = searchsortedfirst(rata2datetime.(time_waves), T_END_ITP)
    time_waves = time_waves[T_START_ITP_idx:T_END_ITP_idx]

    # swh

    swh_waves = ncread(waves_raw, "swh")[:, :, T_START_ITP_idx:T_END_ITP_idx]
    scale_factor = ncgetatt(waves_raw, "swh", "scale_factor")
    units_factor = ncgetatt(waves_raw, "swh", "units") |> x -> uconvert(UNITS_OUT_HEIGHT, 1.0 * NetCDF_units[x]).val
    fill_value = ncgetatt(waves_raw, "swh", "_FillValue")

    swh_waves[swh_waves .== fill_value] .= 0.0
    swh_waves = swh_waves .|> x -> scale_factor * units_factor * x 

    # ust

    ust_waves = ncread(waves_raw, "ust")[:, :, T_START_ITP_idx:T_END_ITP_idx]
    scale_factor = ncgetatt(waves_raw, "ust", "scale_factor")
    units_factor = ncgetatt(waves_raw, "ust", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
    fill_value = ncgetatt(waves_raw, "ust", "_FillValue")

    ust_waves[ust_waves .== fill_value] .= 0.0
    ust_waves = ust_waves .|> x -> scale_factor * units_factor * x 

    # vst

    vst_waves = ncread(waves_raw, "vst")[:, :, T_START_ITP_idx:T_END_ITP_idx]
    scale_factor = ncgetatt(waves_raw, "vst", "scale_factor")
    units_factor = ncgetatt(waves_raw, "vst", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
    fill_value = ncgetatt(waves_raw, "vst", "_FillValue")

    vst_waves[vst_waves .== fill_value] .= 0.0
    vst_waves = vst_waves .|> x -> scale_factor * units_factor * x 

    ####################################################################

    windout = joinpath(@__DIR__, "wind-2018.mat")
    tempout = joinpath(@__DIR__, "temp-2018.mat")
    no3out = joinpath(@__DIR__, "no3-2018.mat")
    waterout = joinpath(@__DIR__, "water-2018.mat")
    wavesout = joinpath(@__DIR__, "waves-2018.mat")

    rm(windout, force = true)
    rm(tempout, force = true)
    rm(no3out, force = true)
    rm(waterout, force = true)
    rm(wavesout, force = true)

    matwrite(windout, Dict("lon" => lon_wind, "lat" => lat_wind, "t" => time_wind, "u" => u_wind, "v" => v_wind))
    matwrite(tempout, Dict("lon" => lon_temp, "lat" => lat_temp, "t" => time_temp, "temp" => temp))
    matwrite(no3out, Dict("lon" => lon_nutr, "lat" => lat_nutr, "t" => time_nutr, "no3" => no3))
    matwrite(waterout, Dict("lon" => lon_water, "lat" => lat_water, "t" => time_water, "u" => u_water, "v" => v_water))
    matwrite(wavesout, Dict("lon" => lon_waves, "lat" => lat_waves, "t" => time_waves, "swh" => swh_waves, "u" => ust_waves, "v" => vst_waves))

    return nothing
end
