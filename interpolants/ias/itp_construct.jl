using JLD2

include(joinpath(@__DIR__, "..", "interpolant-constructors.jl"))

#############################################################################

"""
    const wind_file_ias = "wind-2021-ias.mat"

This is wind velocity data from 2021 in the north Atlantic. The dataset contains 5 keys:

"Lon" [399 x 1] [-99.75:0.25:-0.25 degrees]: The longitude (East/West) at which the measurement is taken.
"Lat" [360 x 1] [90.0:-0.25:0.25 degrees]: The latitude (North/South) at which the measurement is taken.
"t" [364 x 1] [738157.5:1.0:738520.5 Rata Die days}]: The time at which the measurement is taken. See `Dates.rata2datetime`.
"u" [360 x 399 x 364] [lat, lon, time] [km/day]: The x component (East/West) of the wind velocity.
"v" [360 x 399 x 364] [lat, lon, time] [km/day]: The y component (North/South) of the wind velocity.
"""
const wind_file_ias = joinpath(@__DIR__, "wind-2021-ias.mat")

"""
    const water_file_ias = "water-2021-ias.mat"

This is water velocity data from 2021 in the north Atlantic. The dataset contains 5 keys:

"lon" [200 x 1] [-99.875:0.25:-50.125 degrees]: The longitude (East/West) at which the measurement is taken.
"lat" [120 x 1] [5.125:0.25:34.875 degrees]: The latitude (North/South) at which the measurement is taken.
"t" [365 x 1] [738157.0:1.0:738521.0 Rata Die days}]: The time at which the measurement is taken. See `Dates.rata2datetime`.
"u" [120 x 200 x 365] [lat, lon, time] [km/day]: The x component (East/West) of the wind velocity.
"v" [120 x 200 x 365] [lat, lon, time] [km/day]: The y component (North/South) of the wind velocity.
"""
const water_file_ias = joinpath(@__DIR__, "water-2021-ias.mat")

"""
    const ref_ias
    
The default [`EquirectangularReference`](@ref) for the tropical Atlantic.

### Fields
-`lon0`: -75.0 degrees
-`lat0`: 10.0 degrees
-`R`: 6731 km (default)
"""
const ref_ias = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

##############################################################################
##############################################################################
##############################################################################

@info "Constructing wind interpolant."

outfile = joinpath(@__DIR__, "wind_itp.jld2") 
rm(outfile, force = true)

wind_itp = VectorField2DGridSPH(wind_file_ias, lon_alias = "Lon", lat_alias = "Lat", lon_lat_time_order = [2, 1, 3])
wind_itp = VectorField2DInterpolantEQR(wind_itp, ref_ias)
jldsave(outfile, wind_itp = wind_itp)

@info "Wind interpolant written to $(outfile)."

##############################################################################

@info "Constructing water interpolant."

outfile = joinpath(@__DIR__, "water_itp.jld2") 
rm(outfile, force = true)

water_itp = VectorField2DGridSPH(water_file_ias, lon_lat_time_order = [2, 1, 3])
water_itp = VectorField2DInterpolantEQR(water_itp, ref_ias)
jldsave(outfile, water_itp = water_itp)

@info "Water interpolant written to $(outfile)."