using GeoDatasets
using MAT
using JLD2

include(joinpath(@__DIR__, "..", "interpolant-constructors.jl"))

#############################################################################

"""
    const ref_land
    
The default [`EquirectangularReference`](@ref) for the tropical Atlantic.

### Fields
-`lon0`: -75.0 degrees
-`lat0`: 10.0 degrees
-`R`: 6731 km (default)
"""
const ref_land = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

#############################################################################

@info "Constructing land interpolant."

# 0 is ocean, 1 is land and 2 is lake
lon, lat, data = GeoDatasets.landseamask(resolution = 'c', grid = 5)
data[data .== 2] .= 1 # lake is not ocean, so it's land

land_itp = StaticField2DInterpolantEQR(lon, lat, data, ref_land)

outfile = joinpath(@__DIR__, "land_itp.jld2")
rm(outfile, force = true)
jldsave(outfile, land_itp = land_itp)

@info "Land interpolant written to $(outfile)."





