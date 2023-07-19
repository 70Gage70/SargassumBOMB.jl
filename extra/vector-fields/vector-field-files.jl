using JLD2

include("vector-field-methods.jl")

#######################################

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