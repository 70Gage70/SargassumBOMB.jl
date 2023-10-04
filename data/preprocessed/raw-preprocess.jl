using NetCDF
using MAT
using Dates
using Unitful # ms = "m/s", val = 2.0 * uparse(ms)

####################################################################

raw_path = joinpath(@__DIR__, "..", "raw")
wind_raw = joinpath(raw_path, "wind.nc")
water_temp_raw = joinpath(raw_path, "currents-temperature.nc")
nutr_raw = joinpath(raw_path, "nutrients.nc")

NetCDF_units = Dict(
    "m s**-1" => u"m/s",
    "m s-1" => u"m/s", 
    "degrees_C" => u"°C",
    "mmol m-3" => u"mmol/m^3"
)
const UNITS_OUT_SPEED = u"km/d"
const UNITS_OUT_TEMP = u"°C"
const UNITS_OUT_NUTR = u"mmol/m^3"

const t_start = DateTime(2018, 1, 1)
const t_end = DateTime(2018, 12, 31)

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
t_start_idx = searchsortedfirst(rata2datetime.(time_wind), t_start)
t_end_idx = searchsortedfirst(rata2datetime.(time_wind), t_end)
time_wind = time_wind[t_start_idx:t_end_idx]

# u

u_wind = ncread(wind_raw, "u10")[:, :, t_start_idx:t_end_idx]
scale_factor = ncgetatt(wind_raw, "u10", "scale_factor")
units_factor = ncgetatt(wind_raw, "u10", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(wind_raw, "u10", "_FillValue")

u_wind[u_wind .== fill_value] .= 0.0
u_wind = u_wind .|> x -> scale_factor * units_factor * x 

# v

v_wind = ncread(wind_raw, "v10")[:, :, t_start_idx:t_end_idx]
scale_factor = ncgetatt(wind_raw, "v10", "scale_factor")
units_factor = ncgetatt(wind_raw, "v10", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(wind_raw, "v10", "_FillValue")

v_wind[v_wind .== fill_value] .= 0.0
v_wind = v_wind .|> x -> scale_factor * units_factor * x 

###########################
########################### WATER
###########################

# variables: longitude, latitude, depth, time, uo_glor(longitude, latitude, depth time), vo_glor(longitude, latitude, depth, time)
# longitude is in degrees E/W [-180, 180]
# latitude is in degrees N/S [-90 - 90]
# depth is in m
# time is in days since 1950-01-01 00:00:00
# uo_glor and vo_glor are in m/s

lon_water = ncread(water_temp_raw, "longitude") .|> Float64 
lat_water = ncread(water_temp_raw, "latitude") .|> Float64

tref = DateTime(1950, 1, 1, 0, 0, 0)
time_water = ncread(water_temp_raw, "time") .|> (x -> tref + Day(x)) .|> datetime2rata
t_start_idx = searchsortedfirst(rata2datetime.(time_water), t_start)
t_end_idx = searchsortedfirst(rata2datetime.(time_water), t_end)
time_water = time_water[t_start_idx:t_end_idx]

# u

u_water = ncread(water_temp_raw, "uo_glor")[:, :, 1, t_start_idx:t_end_idx]
scale_factor = ncgetatt(water_temp_raw, "uo_glor", "scale_factor")
units_factor = ncgetatt(water_temp_raw, "uo_glor", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(water_temp_raw, "uo_glor", "_FillValue")

u_water[u_water .== fill_value] .= 0.0
u_water = u_water .|> x -> scale_factor * units_factor * x 

# v

v_water = ncread(water_temp_raw, "vo_glor")[:, :, 1, t_start_idx:t_end_idx]
scale_factor = ncgetatt(water_temp_raw, "vo_glor", "scale_factor")
units_factor = ncgetatt(water_temp_raw, "vo_glor", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(water_temp_raw, "vo_glor", "_FillValue")

v_water[v_water .== fill_value] .= 0.0
v_water = v_water .|> x -> scale_factor * units_factor * x 

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
t_start_idx = searchsortedfirst(rata2datetime.(time_temp), t_start)
t_end_idx = searchsortedfirst(rata2datetime.(time_temp), t_end)
time_temp = time_temp[t_start_idx:t_end_idx]

# temp

temp = ncread(water_temp_raw, "thetao_glor")[:, :, 1, t_start_idx:t_end_idx]
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
t_start_idx = searchsortedfirst(rata2datetime.(time_nutr), t_start)
t_end_idx = searchsortedfirst(rata2datetime.(time_nutr), t_end)
time_nutr = time_nutr[t_start_idx:t_end_idx]

# no3

no3 = ncread(nutr_raw, "no3")[:, :, 1, t_start_idx:t_end_idx]
scale_factor = 1.0
units_factor = ncgetatt(nutr_raw, "no3", "units") |> x -> uconvert(UNITS_OUT_NUTR, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(nutr_raw, "no3", "_FillValue")

no3[no3 .== fill_value] .= 0.0
no3 = no3 .|> x -> scale_factor * units_factor * x

####################################################################

windout = joinpath(@__DIR__, "wind-2018.mat")
waterout = joinpath(@__DIR__, "water-2018.mat")
tempout = joinpath(@__DIR__, "temp-2018.mat")
no3out = joinpath(@__DIR__, "no3-2018.mat")

rm(windout, force = true)
rm(waterout, force = true)
rm(tempout, force = true)
rm(no3out, force = true)

matwrite(windout, Dict("lon" => lon_wind, "lat" => lat_wind, "t" => time_wind, "u" => u_wind, "v" => v_wind))
matwrite(waterout, Dict("lon" => lon_water, "lat" => lat_water, "t" => time_water, "u" => u_water, "v" => v_water))
matwrite(tempout, Dict("lon" => lon_temp, "lat" => lat_temp, "t" => time_temp, "temp" => temp))
matwrite(no3out, Dict("lon" => lon_nutr, "lat" => lat_nutr, "t" => time_nutr, "no3" => no3))

include(joinpath(@__DIR__), "rick-preprocess.jl")