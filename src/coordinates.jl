"""
    EquirectangularReference{T}

A container for the reference longitude and latitude of an equirectangular projection.

### Fields

- `lon0`: The standard longitude degrees (East/West).
- `lat0`: The standard latitude degrees (North/South).
"""
struct EquirectangularReference{T<:Real}
    lon0::T
    lat0::T
end


function EquirectangularReference(; lon0::T, lat0::T) where {T<:Real}
    @assert -180.0 <= lon0 <= 180.0 "The longitude must be between -180 degrees and 180 degrees."
    @assert -90 <= lat0 <= 90 "The latitude must be between -90 degrees and 90 degrees."

    return EquirectangularReference(lon0, lat0)
end


"""
    sph2xy(lon, lat, ref; R = 6371)

Compute planar coordinates `[x, y]` [km] from spherical coordinates `(lon, lat)` [deg] with reference `ref::EquirectangularReference`.

### Arguments

`lon`: Longitude in degrees (East/West).
`lat`: Latitude in degrees (North/South).
`ref`: An [`EquirectangularReference`](@ref).

### Optional Arguments

`R`: The radius of the Earth; x and y are returned in the units of R. Default `6371 km`.
"""
function sph2xy(lon, lat, ref; R = 6371)
    @assert -180.0 <= lon <= 180.0 "The longitude must be between -180 degrees and 180 degrees."
    @assert -90 <= lat <= 90 "The latitude must be between -90 degrees and 90 degrees."

    lon0, lat0 = (ref.lon0, ref.lat0)
    deg2rad = π/180

    x = R*(lon - lon0)*deg2rad*cos(lat0*deg2rad)
    y = R*(lat - lat0)*deg2rad

    return [x, y]
end

"""
    xy2sph(x, y, ref; R = 6371)

Compute spherical coordinates `[lon, lat]` [deg] from rectilinear coordinates `(x, y)` [km] from with references `ref::EquirectangularReference`.

### Arguments

`x`: The x Cartesian coordinate in km (East/West).
`y`: The y Cartesian coordinate in km (North/South).
`ref`: An [`EquirectangularReference`](@ref).

### Optional Arguments 

`R`: The radius of the Earth; should be in the same units as x and y. Default 6371 km.
"""
function xy2sph(x, y, ref; R = 6371)
    lon0, lat0 = (ref.lon0, ref.lat0)
    deg2rad = π/180
    rad2deg = 1/deg2rad

    lon = lon0 + rad2deg*x/(R*cos(lat0*deg2rad))
    lat = lat0 + rad2deg*y/R 

    return [lon, lat]
end