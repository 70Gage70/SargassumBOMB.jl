using JLD2

include(joinpath(@__DIR__, "..", "itp-core.jl"))
include(joinpath(@__DIR__, "..", "itp-derivatives.jl"))

#############################################################################

waves_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "waves-2018.mat")

##############################################################################

@info "Constructing waves interpolant."

outfile = joinpath(@__DIR__, "waves_itp.jld2")
rm(outfile, force = true)

itp = GriddedField(waves_file, ["lon", "lat", "t"], ["swh"], 
    time_index = 3, 
    time2datetime = rata2datetime_minute, 
    NaN_replacement = 0.0, 
    var_units = ["deg E/W", "deg N/S", "days"], 
    field_units = ["m"], 
    ref = ref_itp)
itp = itp |> sph2xy |> interpolate

jldsave(outfile, waves_itp = itp)

@info "Waves interpolant written to $(outfile)."

##############################################################################

@info "Constructing Stokes interpolant."

outfile = joinpath(@__DIR__, "stokes_itp.jld2")
rm(outfile, force = true)

itp = GriddedField(waves_file, ["lon", "lat", "t"], ["u", "v"], 
    time_index = 3, 
    time2datetime = rata2datetime_minute, 
    NaN_replacement = 0.0, 
    var_units = ["deg E/W", "deg N/S", "days"], 
    field_units = ["km/d", "km/d"], 
    ref = ref_itp)
itp = itp |> sph2xy |> interpolate
itp = add_derivatives(itp)

jldsave(outfile, waves_itp = itp)

@info "Stokes interpolant written to $(outfile)."

##############################################################################
