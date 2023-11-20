"""
    construct_waves_itp()
"""
function construct_waves_itp()
    @info "Constructing waves interpolant."

    waves_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "waves-2018.mat")
    outfile = joinpath(@__DIR__, "WAVES_ITP.jld2")
    rm(outfile, force = true)

    itp = GriddedField(waves_file, ["lon", "lat", "t"], ["swh"], 
        time_index = 3, 
        time2datetime = rata2datetime_minute, 
        NaN_replacement = 0.0, 
        var_units = ["deg E/W", "deg N/S", "days"], 
        field_units = ["m"], 
        ref = EQR_DEFAULT)
    itp = itp |> sph2xy |> interpolate

    jldsave(outfile, WAVES_ITP = itp)

    @info "Waves interpolant written to $(outfile)."

    return nothing
end

"""
    construct_stokes_itp()
"""
function construct_stokes_itp()
    @info "Constructing Stokes interpolant."

    waves_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "waves-2018.mat")
    outfile = joinpath(@__DIR__, "STOKES_ITP.jld2")
    rm(outfile, force = true)

    itp = GriddedField(waves_file, ["lon", "lat", "t"], ["u", "v"], 
        time_index = 3, 
        time2datetime = rata2datetime_minute, 
        NaN_replacement = 0.0, 
        var_units = ["deg E/W", "deg N/S", "days"], 
        field_units = ["km/d", "km/d"], 
        ref = EQR_DEFAULT)
    itp = itp |> sph2xy |> interpolate
    itp = add_derivatives(itp)

    jldsave(outfile, STOKES_ITP = itp)

    @info "Stokes interpolant written to $(outfile)."

    return nothing
end
