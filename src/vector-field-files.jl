using JLD2

include("vector-field-methods.jl")

#######################################

"""
    wind_file_default = "viento_2021.mat"

This is wind velocity data from 2021 in the north Atlantic. The dataset contains 5 keys:

"Lon" [399 x 1] [-99.75:0.25:-0.25 degrees]: The longitude (East/West) at which the measurement is taken.
"Lat" [360 x 1] [90.0:-0.25:0.25 degrees]: The latitude (North/South) at which the measurement is taken.
"t" [364 x 1] [738157.5:1.0:738520.5 Rata Die days}]: The time at which the measurement is taken. See `Dates.rata2datetime`.
"u" [360 x 399 x 364] [lat, lon, time] [km/day]: The x component (East/West) of the wind velocity.
"v" [360 x 399 x 364] [lat, lon, time] [km/day]: The y component (North/South) of the wind velocity.
"""
const wind_file_default = joinpath(@__DIR__, "..", "MATLAB", "viento_2021.mat")


"""
    water_file_default = "merged-2021-IAS.mat"

This is water velocity data from 2021 in the north Atlantic. The dataset contains 5 keys:

"lon" [200 x 1] [-99.875:0.25:-50.125 degrees]: The longitude (East/West) at which the measurement is taken.
"lat" [120 x 1] [5.125:0.25:34.875 degrees]: The latitude (North/South) at which the measurement is taken.
"t" [365 x 1] [738157.0:1.0:738521.0 Rata Die days}]: The time at which the measurement is taken. See `Dates.rata2datetime`.
"u" [120 x 200 x 365] [lat, lon, time] [km/day]: The x component (East/West) of the wind velocity.
"v" [120 x 200 x 365] [lat, lon, time] [km/day]: The y component (North/South) of the wind velocity.
"""
const water_file_default = joinpath(@__DIR__, "..", "MATLAB", "merged-2021-IAS.mat")

"""
    construct_wind_itp_EQR(infile, ref; outfile)

Build an interpolant for the wind data in `infile` on the equirectangular projection defined by `ref`. 
    
Save the result to a `.jld` file with the name `outfile`.
"""
function construct_wind_itp_EQR(infile::String, ref::EquirectangularReference; outfile::String = "wind_itp.jld")
    @info "Constructing wind interpolant."
    wind_itp = VectorField2DGridSPH(infile, lon_alias = "Lon", lat_alias = "Lat", lon_lat_time_order = [2, 1, 3])
    wind_itp = VectorField2DInterpolantEQR(wind_itp, ref)
    @save outfile wind_itp
    @info "Wind interpolant written to $(outfile)."

    return nothing
end

"""
    construct_water_itp_EQR(infile, ref; outfile)

Build an interpolant for the wind data in `infile` on the equirectangular projection defined by `ref`. 
    
Save the result to a `.jld` file with the name `outfile`.
"""
function construct_water_itp_EQR(infile::String,ref::EquirectangularReference; outfile::String = "water_itp.jld")
    @info "Constructing water interpolant."
    water_itp = VectorField2DGridSPH(infile, lon_lat_time_order = [2, 1, 3])
    water_itp = VectorField2DInterpolantEQR(water_itp, ref)
    @save outfile water_itp
    @info "Water interpolant written to $(outfile)."

    return nothing
end

##############

function construct_itp_EQR(
    water_infile::String, 
    wind_infile::String, 
    ref::EquirectangularReference, 
    clump_parameters::ClumpParameters; 
    outfile::String = "itp.jld2")

    @info "Constructing wind interpolant."
    wind_itp = VectorField2DGridSPH(wind_infile, lon_alias = "Lon", lat_alias = "Lat", lon_lat_time_order = [2, 1, 3])
    wind_itp = VectorField2DInterpolantEQR(wind_itp, ref)
    
    @info "Constructing water interpolant."
    water_itp = VectorField2DGridSPH(water_infile, lon_lat_time_order = [2, 1, 3])
    water_itp = VectorField2DInterpolantEQR(water_itp, ref)
    
    @info "Constructing DDt."
    MDX = MaterialDerivativeX(water_itp)
    MDY = MaterialDerivativeY(water_itp)
    V = Vorticity(water_itp)

    @info "Constructing U."
    UX = WindWaterAlphaX(wind_itp, water_itp, clump_parameters)    
    UY = WindWaterAlphaY(wind_itp, water_itp, clump_parameters) 
    DUX = DDtWindWaterAlphaX(wind_itp, water_itp, clump_parameters)
    DUY = DDtWindWaterAlphaY(wind_itp, water_itp, clump_parameters)

    jldsave(outfile;
        v_x = water_itp.u,
        v_y =  water_itp.v,
        Dv_xDt = MDX,
        Dv_yDt = MDY,
        u_x = UX,
        u_y = UY,
        Du_xDt = DUX, 
        Du_yDt = DUY,
        ω = V
    )

    @info "Interpolants written to $(outfile)."

    return nothing

end