"""
    construct_temp_itp()
"""
function construct_temp_itp()
    @info "Constructing temperature interpolant."

    temp_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "temp-2018.mat")
    outfile = joinpath(@__DIR__, "TEMP_ITP.jld2")
    rm(outfile, force = true)

    itp = GriddedField(temp_file, ["lon", "lat", "t"], ["temp"], 
        time_index = 3, 
        time2datetime = rata2datetime_minute, 
        NaN_replacement = 0.0, 
        var_units = ["deg E/W", "deg N/S", "days"], 
        field_units = ["° C"], 
        ref = EQR.x)
    itp = itp |> sph2xy |> interpolate

    jldsave(outfile, TEMP_ITP = itp)

    @info "Temperature interpolant written to $(outfile)."

    return nothing
end

"""
    construct_no3_itp()
"""
function construct_no3_itp()
    @info "Constructing NO3 interpolant."

    no3_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "no3-2018.mat")
    outfile = joinpath(@__DIR__, "NO3_ITP.jld2")
    rm(outfile, force = true)

    itp = GriddedField(no3_file, ["lon", "lat", "t"], ["no3"], 
        time_index = 3, 
        time2datetime = rata2datetime_minute, 
        NaN_replacement = 0.0, 
        var_units = ["deg E/W", "deg N/S", "days"], 
        field_units = ["mmol/m^3"], 
        ref = EQR.x)
    itp = itp |> sph2xy |> interpolate

    jldsave(outfile, NO3_ITP = itp)

    @info "NO3 interpolant written to $(outfile)."

    return nothing
end
