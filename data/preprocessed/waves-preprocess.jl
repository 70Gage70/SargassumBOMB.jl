using NetCDF
using MAT
using Dates
using Unitful # ms = "m/s", val = 2.0 * uparse(ms)

####################################################################

raw_path = joinpath(@__DIR__, "..", "raw")
waves_raw = joinpath(raw_path, "waves.nc")

NetCDF_units = Dict(
    "m" => u"m",
    "m s**-1" => u"m/s",
    "m s-1" => u"m/s", 
    "degrees_C" => u"°C",
    "mmol m-3" => u"mmol/m^3"
)
const UNITS_OUT_HEIGHT = u"m"
const UNITS_OUT_SPEED = u"km/d"
const UNITS_OUT_TEMP = u"°C"
const UNITS_OUT_NUTR = u"mmol/m^3"

const t_start = DateTime(2018, 1, 1)
const t_end = DateTime(2018, 12, 31)

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
t_start_idx = searchsortedfirst(rata2datetime.(time_waves), t_start)
t_end_idx = searchsortedfirst(rata2datetime.(time_waves), t_end)
time_waves = time_waves[t_start_idx:t_end_idx]

# swh

swh_waves = ncread(waves_raw, "swh")[:, :, t_start_idx:t_end_idx]
scale_factor = ncgetatt(waves_raw, "swh", "scale_factor")
units_factor = ncgetatt(waves_raw, "swh", "units") |> x -> uconvert(UNITS_OUT_HEIGHT, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(waves_raw, "swh", "_FillValue")

swh_waves[swh_waves .== fill_value] .= 0.0
swh_waves = swh_waves .|> x -> scale_factor * units_factor * x 

# ust

ust_waves = ncread(waves_raw, "ust")[:, :, t_start_idx:t_end_idx]
scale_factor = ncgetatt(waves_raw, "ust", "scale_factor")
units_factor = ncgetatt(waves_raw, "ust", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(waves_raw, "ust", "_FillValue")

ust_waves[ust_waves .== fill_value] .= 0.0
ust_waves = ust_waves .|> x -> scale_factor * units_factor * x 

# vst

vst_waves = ncread(waves_raw, "vst")[:, :, t_start_idx:t_end_idx]
scale_factor = ncgetatt(waves_raw, "vst", "scale_factor")
units_factor = ncgetatt(waves_raw, "vst", "units") |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = ncgetatt(waves_raw, "vst", "_FillValue")

vst_waves[vst_waves .== fill_value] .= 0.0
vst_waves = vst_waves .|> x -> scale_factor * units_factor * x 

####################################################################

wavesout = joinpath(@__DIR__, "waves-2018.mat")

rm(wavesout, force = true)

matwrite(wavesout, Dict("lon" => lon_waves, "lat" => lat_waves, "t" => time_waves, "swh" => swh_waves, "u" => ust_waves, "v" => vst_waves))