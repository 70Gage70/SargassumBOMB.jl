using NetCDF
using MAT
using Dates
using Unitful # ms = "m/s", val = 2.0 * uparse(ms)

####################################################################

raw_path = joinpath(@__DIR__, "..", "raw-data")
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
const UNITS_OUT_NUTR = u"Mmol/km^3"

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
dec1_2017 = searchsortedfirst(rata2datetime.(time_wind), DateTime(2017, 12, 1))
jan31_2019 = searchsortedfirst(rata2datetime.(time_wind), DateTime(2019, 1, 31))
time_wind = time_wind[dec1_2017:jan31_2019]

# u

u_wind = ncread(wind_raw, "u10")[:,:,dec1_2017:jan31_2019]
scale_factor = ncgetatt(wind_raw, "u10", "scale_factor")
units_factor = ncgetatt(wind_raw, "u10", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(wind_raw, "u10", "_FillValue")

u_wind[u_wind .== fill_value] .= 0.0
u_wind = u_wind .|> x -> scale_factor * units_factor * x 

# v

v_wind = ncread(wind_raw, "v10")[:,:,dec1_2017:jan31_2019]
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

# u

u_water = ncread(water_temp_raw, "uo_glor")[:, :, 1, :]
scale_factor = ncgetatt(water_temp_raw, "uo_glor", "scale_factor")
units_factor = ncgetatt(water_temp_raw, "uo_glor", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(water_temp_raw, "uo_glor", "_FillValue")

u_water[u_water .== fill_value] .= 0.0
u_water = u_water .|> x -> scale_factor * units_factor * x 

# v

v_water = ncread(water_temp_raw, "vo_glor")[:, :, 1, :]
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

# temp

temp = ncread(water_temp_raw, "thetao_glor")[:, :, 1, :]
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

# no3

no3 = ncread(nutr_raw, "no3")[:, :, 1, :]
scale_factor = 1.0
units_factor = ncgetatt(nutr_raw, "no3", "units") |> x -> uconvert(UNITS_OUT_NUTR, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(nutr_raw, "no3", "_FillValue")

no3[no3 .== fill_value] .= 0.0
no3 = no3 .|> x -> scale_factor * units_factor * x

####################################################################

rm("wind-2018.mat", force = true)
rm("water-2018.mat", force = true)
rm("temp-2018.mat", force = true)
rm("no3-2018.mat", force = true)

matwrite("wind-2018.mat", Dict("lon" => lon_wind, "lat" => lat_wind, "t" => time_wind, "u" => u_wind, "v" => v_wind))
matwrite("water-2018.mat", Dict("lon" => lon_water, "lat" => lat_water, "t" => time_water, "u" => u_water, "v" => v_water))
matwrite("temp-2018.mat", Dict("lon" => lon_temp, "lat" => lat_temp, "t" => time_temp, "u" => temp))
matwrite("no3-2018.mat", Dict("lon" => lon_nutr, "lat" => lat_nutr, "t" => time_nutr, "u" => no3))