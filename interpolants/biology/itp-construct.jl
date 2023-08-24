using JLD2

include(joinpath(@__DIR__, "..", "itp-core.jl"))

#############################################################################

temp_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "temp-2018.mat")
no3_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "no3-2018.mat")

##############################################################################

@info "Constructing temperature interpolant."

outfile = joinpath(@__DIR__, "temp_itp.jld2")
rm(outfile, force = true)

itp = GriddedField(temp_file, ["lon", "lat", "t"], ["temp"], 
    time_index = 3, 
    time2datetime = rata2datetime_minute, 
    NaN_replacement = 0.0, 
    var_units = ["deg E/W", "deg N/S", "days"], 
    field_units = ["° C"], 
    ref = ref_itp)
itp = itp |> sph2xy |> interpolate

jldsave(outfile, temp_itp = itp)

@info "Temperature interpolant written to $(outfile)."

##############################################################################

@info "Constructing NO3 interpolant."

outfile = joinpath(@__DIR__, "no3_itp.jld2")
rm(outfile, force = true)

itp = GriddedField(no3_file, ["lon", "lat", "t"], ["no3"], 
    time_index = 3, 
    time2datetime = rata2datetime_minute, 
    NaN_replacement = 0.0, 
    var_units = ["deg E/W", "deg N/S", "days"], 
    field_units = ["mmol/m^3"], 
    ref = ref_itp)
itp = itp |> sph2xy |> interpolate

jldsave(outfile, no3_itp = itp)

@info "NO3 interpolant written to $(outfile)."
