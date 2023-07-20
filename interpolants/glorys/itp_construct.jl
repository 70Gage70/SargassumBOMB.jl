using JLD2

include(joinpath(@__DIR__, "..", "interpolant-constructors.jl"))

#############################################################################

"""
    const wind_file_glor = "wind-2018-glor.mat"
"""
const wind_file_glor = joinpath(@__DIR__, "wind-2018-glor.mat")

"""
    const water_file_glor = "water-2018-glor.mat"
"""
const water_file_glor = joinpath(@__DIR__, "water-2018-glor.mat")

"""
    const temp_file_glor = "temp-2018-glor.mat"
"""
const temp_file_glor = joinpath(@__DIR__, "temp-2018-glor.mat")

"""
    const no3_file_glor = "no3-2018-glor.mat"
"""
const no3_file_glor = joinpath(@__DIR__, "no3-2018-glor.mat")

"""
    const ref_glor
    
The default [`EquirectangularReference`](@ref) for the tropical Atlantic.

### Fields
-`lon0`: -75.0 degrees
-`lat0`: 10.0 degrees
-`R`: 6731 km (default)
"""
const ref_glor = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

##############################################################################
##############################################################################
##############################################################################

@info "Constructing wind interpolant."

outfile = joinpath(@__DIR__, "wind_itp.jld2")
rm(outfile, force = true)

itp = VectorField2DGridSPH(wind_file_glor)
itp = VectorField2DInterpolantEQR(itp, ref_glor)
jldsave(outfile, wind_itp = itp)

@info "Wind interpolant written to $(outfile)."

##############################################################################

@info "Constructing water interpolant."

outfile = joinpath(@__DIR__, "water_itp.jld2")
rm(outfile, force = true)

itp = VectorField2DGridSPH(water_file_glor)
itp = VectorField2DInterpolantEQR(itp, ref_glor)
jldsave(outfile, water_itp = itp)

@info "Water interpolant written to $(outfile)."

##############################################################################

@info "Constructing temperature interpolant."

outfile = joinpath(@__DIR__, "temp_itp.jld2")
rm(outfile, force = true)

itp = ScalarField2DGridSPH(temp_file_glor)
itp = ScalarField2DInterpolantEQR(itp, ref_glor)
jldsave(outfile, temp_itp = itp)

@info "Temperature interpolant written to $(outfile)."

##############################################################################

@info "Constructing NO3 interpolant."

outfile = joinpath(@__DIR__, "no3_itp.jld2")
rm(outfile, force = true)

itp = ScalarField2DGridSPH(no3_file_glor)
itp = ScalarField2DInterpolantEQR(itp, ref_glor)
jldsave(outfile, no3_itp = itp)

@info "NO3 interpolant written to $(outfile)."