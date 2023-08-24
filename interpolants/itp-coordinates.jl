include(joinpath(@__DIR__, "..", "src", "coordinates.jl"))

#####################################################################################

"""
    const ref_itp
    
The default [`EquirectangularReference`](@ref) for the tropical Atlantic.

### Fields
-`lon0`: -75.0 degrees
-`lat0`: 10.0 degrees
-`R`: 6731 km (default)
"""
const ref_itp = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)