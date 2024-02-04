"""
    construct_wind_itp()
"""
function construct_wind_itp()
    @info "Constructing wind interpolant."

    wind_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "wind-2018.mat")
    outfile = joinpath(@__DIR__, "WIND_ITP.jld2")
    rm(outfile, force = true)

    itp = GriddedField(wind_file, ["lon", "lat", "t"], ["u", "v"], 
        time_index = 3, 
        time2datetime = rata2datetime_minute, 
        NaN_replacement = 0.0, 
        var_units = ["deg E/W", "deg N/S", "days"], 
        field_units = ["km/d", "km/d"], 
        ref = EQR_DEFAULT)
    itp = itp |> sph2xy |> interpolate
    itp = add_derivatives(itp)

    jldsave(outfile, WIND_ITP = itp)

    @info "Wind interpolant written to $(outfile)."

    return nothing
end

"""
    construct_water_itp()
"""
function construct_water_itp()
    @info "Constructing water interpolant."

    water_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "water-2018.mat")
    outfile = joinpath(@__DIR__, "WATER_ITP.jld2")
    rm(outfile, force = true)

    itp = GriddedField(water_file, ["lon", "lat", "t"], ["u", "v"], 
        time_index = 3, 
        time2datetime = rata2datetime_minute, 
        NaN_replacement = 0.0, 
        var_units = ["deg E/W", "deg N/S", "days"], 
        field_units = ["km/d", "km/d"], 
        ref = EQR_DEFAULT)
    itp = itp |> sph2xy |> interpolate
    itp = add_derivatives(itp)

    jldsave(outfile, WATER_ITP = itp)

    @info "Water interpolant written to $(outfile)."

    return nothing
end
