using JLD2

include(joinpath(@__DIR__, "..", "interpolant-constructors.jl"))

#############################################################################

"""
    const wind_file = "wind-2018.mat"
"""
const wind_file = joinpath(@__DIR__, "preprocessed", "wind-2018.mat")

"""
    const water_file = "water-2018.mat"
"""
const water_file = joinpath(@__DIR__, "preprocessed", "water-2018.mat")

"""
    const temp_file = "temp-2018.mat"
"""
const temp_file = joinpath(@__DIR__, "preprocessed", "temp-2018.mat")

"""
    const no3_file = "no3-2018.mat"
"""
const no3_file = joinpath(@__DIR__, "preprocessed", "no3-2018.mat")

"""
    const ref_itp
    
The default [`EquirectangularReference`](@ref) for the tropical Atlantic.

### Fields
-`lon0`: -75.0 degrees
-`lat0`: 10.0 degrees
-`R`: 6731 km (default)
"""
const ref_itp = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

##############################################################################
##############################################################################
##############################################################################

@info "Constructing wind interpolant."

outfile = joinpath(@__DIR__, "wind_itp.jld2")
rm(outfile, force = true)

itp = VectorField2DGridSPH(wind_file)
itp = VectorField2DInterpolantEQR(itp, ref_itp)
jldsave(outfile, wind_itp = itp)

@info "Wind interpolant written to $(outfile)."

##############################################################################

@info "Constructing water interpolant."

outfile = joinpath(@__DIR__, "water_itp.jld2")
rm(outfile, force = true)

itp = VectorField2DGridSPH(water_file)
itp = VectorField2DInterpolantEQR(itp, ref_itp)
jldsave(outfile, water_itp = itp)

@info "Water interpolant written to $(outfile)."

##############################################################################

@info "Constructing temperature interpolant."

outfile = joinpath(@__DIR__, "temp_itp.jld2")
rm(outfile, force = true)

itp = ScalarField2DGridSPH(temp_file)
itp = ScalarField2DInterpolantEQR(itp, ref_itp)
jldsave(outfile, temp_itp = itp)

@info "Temperature interpolant written to $(outfile)."

##############################################################################

@info "Constructing NO3 interpolant."

outfile = joinpath(@__DIR__, "no3_itp.jld2")
rm(outfile, force = true)

itp = ScalarField2DGridSPH(no3_file)
itp = ScalarField2DInterpolantEQR(itp, ref_itp)
jldsave(outfile, no3_itp = itp)

@info "NO3 interpolant written to $(outfile)."

# Plotting

# include("../../../CustomMakie.jl/src/geo-methods.jl");
# fig = fig = default_fig();
# ax = geo_axis(fig[1, 1], title = "Wind", limits = (-100, -50, 5, 35));
# wind_itp_test = load("wind_itp.jld2", "wind_itp");