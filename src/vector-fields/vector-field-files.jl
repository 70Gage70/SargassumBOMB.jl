using JLD2

include("vector-field-methods.jl")

#######################################

"""
    const wind_file_default = "viento_2021.mat"

This is wind velocity data from 2021 in the north Atlantic. The dataset contains 5 keys:

"Lon" [399 x 1] [-99.75:0.25:-0.25 degrees]: The longitude (East/West) at which the measurement is taken.
"Lat" [360 x 1] [90.0:-0.25:0.25 degrees]: The latitude (North/South) at which the measurement is taken.
"t" [364 x 1] [738157.5:1.0:738520.5 Rata Die days}]: The time at which the measurement is taken. See `Dates.rata2datetime`.
"u" [360 x 399 x 364] [lat, lon, time] [km/day]: The x component (East/West) of the wind velocity.
"v" [360 x 399 x 364] [lat, lon, time] [km/day]: The y component (North/South) of the wind velocity.
"""
const wind_file_default = joinpath(@__DIR__, "..", "..", "interpolants", "viento_2021.mat")

"""
    const water_file_default = "merged-2021-IAS.mat"

This is water velocity data from 2021 in the north Atlantic. The dataset contains 5 keys:

"lon" [200 x 1] [-99.875:0.25:-50.125 degrees]: The longitude (East/West) at which the measurement is taken.
"lat" [120 x 1] [5.125:0.25:34.875 degrees]: The latitude (North/South) at which the measurement is taken.
"t" [365 x 1] [738157.0:1.0:738521.0 Rata Die days}]: The time at which the measurement is taken. See `Dates.rata2datetime`.
"u" [120 x 200 x 365] [lat, lon, time] [km/day]: The x component (East/West) of the wind velocity.
"v" [120 x 200 x 365] [lat, lon, time] [km/day]: The y component (North/South) of the wind velocity.
"""
const water_file_default = joinpath(@__DIR__, "..", "..", "interpolants", "merged-2021-IAS.mat")

"""
    const ref_default
    
The default [`EquirectangularReference`](@ref) for the tropical Atlantic.
"""
const ref_default = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

"""
    construct_wind_itp_EQR(infile, ref; outfile)

Build an interpolant for the wind data in `infile` on the equirectangular projection defined by `ref`. 
    
Save the result to `outfile` which should be of the form `filename.jld2`.
"""
function construct_wind_itp_EQR(
    infile::String = wind_file_default, 
    ref::EquirectangularReference = ref_default; 
    outfile::String = "wind_itp.jld2")

    @info "Constructing wind interpolant."

    wind_itp = VectorField2DGridSPH(infile, lon_alias = "Lon", lat_alias = "Lat", lon_lat_time_order = [2, 1, 3])
    wind_itp = VectorField2DInterpolantEQR(wind_itp, ref)
    jldsave(outfile, wind_itp = wind_itp)

    @info "Wind interpolant written to $(outfile)."

    return nothing
end

"""
    construct_water_itp_EQR(infile, ref; outfile)

Build an interpolant for the wind data in `infile` on the equirectangular projection defined by `ref`. 
    
    Save the result to `outfile` which should be of the form `filename.jld2`.
"""
function construct_water_itp_EQR(
    infile::String = water_file_default,
    ref::EquirectangularReference = ref_default; 
    outfile::String = "water_itp.jld2")

    @info "Constructing water interpolant."

    water_itp = VectorField2DGridSPH(infile, lon_lat_time_order = [2, 1, 3])
    water_itp = VectorField2DInterpolantEQR(water_itp, ref)
    jldsave(outfile, water_itp = water_itp)

    @info "Water interpolant written to $(outfile)."

    return nothing
end