const UNITS_OUT_HEIGHT = u"m"
const UNITS_OUT_SPEED = u"km/d"
const UNITS_OUT_TEMP = u"°C"
const UNITS_OUT_NUTR = u"mmol/m^3"

const T_START_ITP = DateTime(2018, 1, 1)
const T_END_ITP = DateTime(2018, 12, 31)

"""
    preprocess_all()
"""
function preprocess_all()
    raw_path = joinpath(@__DIR__, "..", "raw")

    wind_raw = joinpath(raw_path, "wind.nc")
    water_temp_raw = joinpath(raw_path, "currents-temperature.nc")
    nutr_raw = joinpath(raw_path, "nutrients.nc")
    rick_raw = matopen(joinpath(raw_path, "merged_2018.mat"))
    waves_raw = joinpath(raw_path, "waves.nc")

    NetCDF_units = Dict(
        "m" => u"m",
        "m s**-1" => u"m/s",
        "m s-1" => u"m/s", 
        "degrees_C" => u"°C",
        "mmol m-3" => u"mmol/m^3"
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
    ########################### RICK'S DATA (WATER)
    ###########################

    # indep vars: "Lat", "Lon", "T" (in that order!)
    # dep vars: "Uek_a", "Uek_bar", "Ug_a", "Ug_bar", "Uslip_d", "Uslip_ud" and similarly for V
    # longitude is in degrees E/W [-180, 180]
    # latitude is in degrees N/S [-90 - 90]
    # time is in Julian days
    # V/U are in m/s

    # to compute the actual water velocities:
    # U = Uek_bar + Uek_a + Ug_bar + Ug_a + Uslip_d (ZONAL - x)
    # V = Vek_bar + Vek_a + Vg_bar + Vg_a + Vslip_d (MERIDIONAL - y)

    lon_rick = read(rick_raw, "Lon") |> vec .|> Float64
    lon_idx = findall(x -> -101 <= x <= -39, lon_rick)
    lon_rick = lon_rick[lon_idx]

    lat_rick = read(rick_raw, "Lat") |> vec .|> Float64
    lat_idx = findall(x -> 0 <= x <= 40, lat_rick)
    lat_rick = lat_rick[lat_idx]

    time_rick = read(rick_raw, "T") |> vec .|> julian2datetime .|> datetime2rata
    T_START_ITP_idx = searchsortedfirst(rata2datetime.(time_rick), T_START_ITP)
    T_END_ITP_idx = searchsortedfirst(rata2datetime.(time_rick), T_END_ITP)
    time_rick = time_rick[T_START_ITP_idx:T_END_ITP_idx]

    # u

    Uek_bar = read(rick_raw, "Uek_bar")[lat_idx, lon_idx] |> x -> permutedims(x)
    Uek_a = read(rick_raw, "Uek_a")[lat_idx, lon_idx, T_START_ITP_idx:T_END_ITP_idx] |> x -> permutedims(x, (2, 1, 3))
    Uek =  zeros(eltype(Uek_a), size(Uek_a))
    for t = 1:size(Uek_a, 3)
        Uek[:,:,t] = Uek_a[:,:,t] + Uek_bar
    end

    Ug_bar = read(rick_raw, "Ug_bar")[lat_idx, lon_idx] |> x -> permutedims(x)
    Ug_a = read(rick_raw, "Ug_a")[lat_idx, lon_idx, T_START_ITP_idx:T_END_ITP_idx] |> x -> permutedims(x, (2, 1, 3))
    Ug =  zeros(eltype(Ug_a), size(Ug_a))
    for t = 1:size(Ug_a, 3)
        Ug[:,:,t] = Ug_a[:,:,t] + Ug_bar
    end

    Uslip_d = read(rick_raw, "Uslip_d")[lat_idx, lon_idx, T_START_ITP_idx:T_END_ITP_idx] |> x -> permutedims(x, (2, 1, 3))

    u_rick = Uek + Ug  + Uslip_d

    scale_factor = 1.0
    units_factor = "m s-1" |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
    fill_value = NaN

    u_rick[u_rick .== fill_value] .= 0.0
    u_rick = u_rick .|> x -> scale_factor * units_factor * x 

    # v

    Vek_bar = read(rick_raw, "Vek_bar")[lat_idx, lon_idx] |> x -> permutedims(x)
    Vek_a = read(rick_raw, "Vek_a")[lat_idx, lon_idx, T_START_ITP_idx:T_END_ITP_idx] |> x -> permutedims(x, (2, 1, 3))
    Vek =  zeros(eltype(Vek_a), size(Vek_a))
    for t = 1:size(Vek_a, 3)
        Vek[:,:,t] = Vek_a[:,:,t] + Vek_bar
    end

    Vg_bar = read(rick_raw, "Vg_bar")[lat_idx, lon_idx] |> x -> permutedims(x)
    Vg_a = read(rick_raw, "Vg_a")[lat_idx, lon_idx, T_START_ITP_idx:T_END_ITP_idx] |> x -> permutedims(x, (2, 1, 3))
    Vg =  zeros(eltype(Vg_a), size(Vg_a))
    for t = 1:size(Vg_a, 3)
        Vg[:,:,t] = Vg_a[:,:,t] + Vg_bar
    end

    Vslip_d = read(rick_raw, "Vslip_d")[lat_idx, lon_idx, T_START_ITP_idx:T_END_ITP_idx] |> x -> permutedims(x, (2, 1, 3))

    v_rick = Vek + Vg  + Vslip_d

    scale_factor = 1.0
    units_factor = "m s-1" |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
    fill_value = NaN

    v_rick[v_rick .== fill_value] .= 0.0
    v_rick = v_rick .|> x -> scale_factor * units_factor * x 

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
    rickout = joinpath(@__DIR__, "rick-2018.mat")
    wavesout = joinpath(@__DIR__, "waves-2018.mat")

    rm(windout, force = true)
    rm(tempout, force = true)
    rm(no3out, force = true)
    rm(rickout, force = true)
    rm(wavesout, force = true)

    matwrite(windout, Dict("lon" => lon_wind, "lat" => lat_wind, "t" => time_wind, "u" => u_wind, "v" => v_wind))
    matwrite(tempout, Dict("lon" => lon_temp, "lat" => lat_temp, "t" => time_temp, "temp" => temp))
    matwrite(no3out, Dict("lon" => lon_nutr, "lat" => lat_nutr, "t" => time_nutr, "no3" => no3))
    matwrite(rickout, Dict("lon" => lon_rick, "lat" => lat_rick, "t" => time_rick, "u" => u_rick, "v" => v_rick))
    matwrite(wavesout, Dict("lon" => lon_waves, "lat" => lat_waves, "t" => time_waves, "swh" => swh_waves, "u" => ust_waves, "v" => vst_waves))

    return nothing
end
