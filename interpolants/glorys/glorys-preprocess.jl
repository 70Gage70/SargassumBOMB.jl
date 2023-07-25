using NetCDF
using MAT
using Dates
using Unitful # ms = "m/s", val = 2.0 * uparse(ms)

####################################################################

# wind_file = joinpath(@__DIR__, "KNMI-GLO-WIND_L3-REP-OBS_METOP-A_ASCAT_25_DES_1688757095019.nc")
wind_file = joinpath(@__DIR__, "cmems_obs-wind_glo_phy_my_l4_P1M_1690305893418.nc")

water_temp_file = joinpath(@__DIR__, "global-reanalysis-phy-001-031-grepv2-daily_1688755986955.nc")

nutr_file = joinpath(@__DIR__, "cmems_mod_glo_bgc_my_0.25_P1D-m_1688755326762.nc")

NOAA_units = Dict(
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

# variables: lon, lat, time, eastward_wind(lon, lat, time), northward_wind(lon, lat, time)
# lon is in degrees East [0 - 360]
# lat is in degrees N/S [-90 - 90]
# time is in seconds since 1990-01-01 00:00:00
# eastward_wind and northward_wind are in m/s

lon_wind = ncread(wind_file, "lon") .|> Float64 # .|> (x -> x - 360.0)
lat_wind = ncread(wind_file, "lat") .|> Float64

# tref = DateTime(1990, 1, 1, 0, 0, 0) |> datetime2unix
# time_wind = ncread(wind_file, "time") .|> (x -> x + tref) .|> unix2datetime .|> datetime2rata
time_wind = [736664 + 30*i for i = 0:13] # for monthy gird

# u

u_wind = ncread(wind_file, "eastward_wind")
scale_factor = ncgetatt(wind_file, "eastward_wind", "scale_factor")
units_factor = ncgetatt(wind_file, "eastward_wind", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NOAA_units[x]).val
fill_value = ncgetatt(wind_file, "eastward_wind", "_FillValue")

u_wind[u_wind .== fill_value] .= 0.0
u_wind = u_wind .|> x -> scale_factor * units_factor * x 

# v

v_wind = ncread(wind_file, "northward_wind")
scale_factor = ncgetatt(wind_file, "northward_wind", "scale_factor")
units_factor = ncgetatt(wind_file, "northward_wind", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NOAA_units[x]).val
fill_value = ncgetatt(wind_file, "northward_wind", "_FillValue")

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

lon_water = ncread(water_temp_file, "longitude") .|> Float64 
lat_water = ncread(water_temp_file, "latitude") .|> Float64

tref = DateTime(1950, 1, 1, 0, 0, 0)
time_water = ncread(water_temp_file, "time") .|> (x -> Day(x) + tref) .|> datetime2rata

# u

u_water = ncread(water_temp_file, "uo_glor")[:, :, 1, :]
scale_factor = ncgetatt(water_temp_file, "uo_glor", "scale_factor")
units_factor = ncgetatt(water_temp_file, "uo_glor", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NOAA_units[x]).val
fill_value = ncgetatt(water_temp_file, "uo_glor", "_FillValue")

u_water[u_water .== fill_value] .= 0.0
u_water = u_water .|> x -> scale_factor * units_factor * x 

# v

v_water = ncread(water_temp_file, "vo_glor")[:, :, 1, :]
scale_factor = ncgetatt(water_temp_file, "vo_glor", "scale_factor")
units_factor = ncgetatt(water_temp_file, "vo_glor", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NOAA_units[x]).val
fill_value = ncgetatt(water_temp_file, "vo_glor", "_FillValue")

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

lon_temp = ncread(water_temp_file, "longitude") .|> Float64 
lat_temp = ncread(water_temp_file, "latitude") .|> Float64

tref = DateTime(1950, 1, 1, 0, 0, 0)
time_temp = ncread(water_temp_file, "time") .|> (x -> Day(x) + tref) .|> datetime2rata

# temp

temp = ncread(water_temp_file, "thetao_glor")[:, :, 1, :]
scale_factor = ncgetatt(water_temp_file, "thetao_glor", "scale_factor")
units_factor = ncgetatt(water_temp_file, "thetao_glor", "units") |> x -> uconvert(UNITS_OUT_TEMP, 1.0 * NOAA_units[x]).val
fill_value = ncgetatt(water_temp_file, "thetao_glor", "_FillValue")

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

lon_nutr = ncread(nutr_file, "longitude") .|> Float64 
lat_nutr = ncread(nutr_file, "latitude") .|> Float64

tref = DateTime(1950, 1, 1, 0, 0, 0)
time_nutr = ncread(nutr_file, "time") .|> (x -> Hour(x) + tref) .|> datetime2rata

# no3

no3 = ncread(nutr_file, "no3")[:, :, 1, :]
scale_factor = 1.0
units_factor = ncgetatt(nutr_file, "no3", "units") |> x -> uconvert(UNITS_OUT_NUTR, 1.0 * NOAA_units[x]).val
fill_value = ncgetatt(nutr_file, "no3", "_FillValue")

no3[no3 .== fill_value] .= 0.0
no3 = no3 .|> x -> scale_factor * units_factor * x

####################################################################

rm("wind-2018-glor.mat", force = true)
rm("water-2018-glor.mat", force = true)
rm("temp-2018-glor.mat", force = true)
rm("no3-2018-glor.mat", force = true)

matwrite("wind-2018-glor.mat", Dict("lon" => lon_wind, "lat" => lat_wind, "t" => time_wind, "u" => u_wind, "v" => v_wind))
matwrite("water-2018-glor.mat", Dict("lon" => lon_water, "lat" => lat_water, "t" => time_water, "u" => u_water, "v" => v_water))
matwrite("temp-2018-glor.mat", Dict("lon" => lon_temp, "lat" => lat_temp, "t" => time_temp, "u" => temp))
matwrite("no3-2018-glor.mat", Dict("lon" => lon_nutr, "lat" => lat_nutr, "t" => time_nutr, "u" => no3))