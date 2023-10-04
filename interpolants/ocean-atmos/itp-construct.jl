using JLD2

include(joinpath(@__DIR__, "..", "itp-core.jl"))
include(joinpath(@__DIR__, "..", "itp-derivatives.jl"))

#############################################################################

wind_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "wind-2018.mat")
water_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "water-2018.mat")
rick_file = joinpath(@__DIR__, "..", "..", "data", "preprocessed", "rick-2018.mat")

##############################################################################

@info "Constructing wind interpolant."

outfile = joinpath(@__DIR__, "wind_itp.jld2")
rm(outfile, force = true)

itp = GriddedField(wind_file, ["lon", "lat", "t"], ["u", "v"], 
    time_index = 3, 
    time2datetime = rata2datetime_minute, 
    NaN_replacement = 0.0, 
    var_units = ["deg E/W", "deg N/S", "days"], 
    field_units = ["km/d", "km/d"], 
    ref = ref_itp)
itp = itp |> sph2xy |> interpolate
itp = add_derivatives(itp)

jldsave(outfile, wind_itp = itp)

@info "Wind interpolant written to $(outfile)."

##############################################################################

@info "Constructing water interpolant."

outfile = joinpath(@__DIR__, "water_itp.jld2")
rm(outfile, force = true)

itp = GriddedField(water_file, ["lon", "lat", "t"], ["u", "v"], 
    time_index = 3, 
    time2datetime = rata2datetime_minute, 
    NaN_replacement = 0.0, 
    var_units = ["deg E/W", "deg N/S", "days"], 
    field_units = ["km/d", "km/d"], 
    ref = ref_itp)
itp = itp |> sph2xy |> interpolate
itp = add_derivatives(itp)

jldsave(outfile, water_itp = itp)

@info "Water interpolant written to $(outfile)."

##############################################################################

@info "Constructing rick interpolant."

outfile = joinpath(@__DIR__, "rick_itp.jld2")
rm(outfile, force = true)

itp = GriddedField(rick_file, ["lon", "lat", "t"], ["u", "v"], 
    time_index = 3, 
    time2datetime = rata2datetime_minute, 
    NaN_replacement = 0.0, 
    var_units = ["deg E/W", "deg N/S", "days"], 
    field_units = ["km/d", "km/d"], 
    ref = ref_itp)
itp = itp |> sph2xy |> interpolate
itp = add_derivatives(itp)

jldsave(outfile, rick_itp = itp)

@info "Rick interpolant written to $(outfile)."
