####################################################################
#################################################################### PREPROCESS
####################################################################

using NetCDF
using MAT
using Dates
using Unitful # ms = "m/s", val = 2.0 * uparse(ms)
using Interpolations

afai_raw = "afai.nc"
ncinfo(afai_raw)

###########################
########################### AFAI
###########################

# variables: longitude, latitude, time, AFAI(longitude, latitude, time)
# longitude is in degrees E/W [-180, 180]
# latitude is in degrees N/S [-90, 90]
# time is in seconds since 1970-01-01T00:00:00
# AFAI is an index measured in reflectance

lon_afai = ncread(afai_raw, "longitude") .|> Float64
lat_afai = ncread(afai_raw, "latitude") .|> Float64

tref = DateTime(1970, 1, 1, 0, 0, 0)
time_afai = ncread(afai_raw, "time") .|> (x -> tref + Second(x)) .|> datetime2rata

# AFAI
# this quantitity has no units or scale factor, and the fill value is NaN


afai = ncread(afai_raw, "AFAI")
fill_value = ncgetatt(afai_raw, "u10", "_FillValue")

afai = afai .|> Float64
afai[isnan.(afai)] .= 0.0

# the data are not completely uniformly gridded, so we we-interpolate on such a grid

gridded_itp = Interpolations.interpolate((lon_afai, lat_afai, time_afai), afai, Gridded(Linear()))
lon_range = range(lon_afai[1], lon_afai[end], step = 0.25) |> collect
lat_range = range(lat_afai[1], lat_afai[end], step = 0.25) |> collect
time_range = range(time_afai[1], time_afai[end], step = 1) |> collect
afai_gridded = [gridded_itp[lon, lat, time] for lon in lon_range, lat in lat_range, time in time_range]

####################################################################

afaiout = "afai-2018.mat"
rm(afaiout, force = true)
matwrite(afaiout, Dict("lon" => lon_range, "lat" => lat_range, "t" => time_range, "afai" => afai_gridded))

####################################################################
#################################################################### Interpolant
####################################################################

using JLD2

include("C:\\Users\\gageb\\Desktop\\Postdoc\\Repositories\\SargassumBOMB.jl\\interpolants\\itp-core.jl")

afai_file = "afai-2018.mat"

@info "Constructing AFAI interpolant."

outfile = joinpath(@__DIR__, "afai_itp.jld2")
rm(outfile, force = true)

itp = GriddedField(afai_file, ["lon", "lat", "t"], ["afai"], 
    time_index = 3, 
    time2datetime = rata2datetime_minute, 
    NaN_replacement = 0.0, 
    var_units = ["deg E/W", "deg N/S", "days"], 
    field_units = ["reflectance"], 
    ref = ref_itp)
itp = itp |> sph2xy |> interpolate

jldsave(outfile, afai_itp = itp)

@info "AFAI interpolant written to $(outfile)."

####################################################################
#################################################################### Plotting
####################################################################

afai_itp = itp

fig = default_fig()
ax = geo_axis(fig[1, 1])

lims = ax.limits.val
n_points = 100
t0 = 121.0 # april 1

xs = range(start = lims[1], stop = lims[2], length = n_points)
ys = range(start = lims[3], stop = lims[4], length = n_points)
zs = [afai_itp.fields[:afai](sph2xy(x, y, afai_itp.ref)..., t) for x in xs, y in ys]
# zs = [sum(afai_itp.fields[:afai](sph2xy(x, y, afai_itp.ref)..., t) for t = t0:t0+31)  for x in xs, y in ys]

heatmap!(ax, xs, ys, zs)

land!(ax)

fig


