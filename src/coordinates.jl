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

"""
    EquirectangularReference(; lon0, lat0)

Construct an `EquirectangularReference` with longitude `lon0` in degrees (East/West) and latitude `lat0` in degrees (North/South).
"""
function EquirectangularReference(; lon0::Real, lat0::Real)
    @assert -180.0 <= lon0 <= 180.0 "The longitude must be between -180 degrees and 180 degrees."
    @assert -90 <= lat0 <= 90 "The latitude must be between -90 degrees and 90 degrees."

    lon0, lat0 = promote(float(lon0), float(lat0))

    return EquirectangularReference(lon0, lat0)
end

"""
    sph2xy(lon, lat, ref; R = 6371)

Compute planar coordinates `[x, y]` from spherical coordinates `(lon, lat)` [deg] with reference `ref::EquirectangularReference`.

The units of `x` and `y` are controlled by the optional argument `R`, the radius of the Earth. The default is `R = 6371 km.`

### Arguments

`lon`: Longitude in degrees (East/West).
`lat`: Latitude in degrees (North/South).
`ref`: An [`EquirectangularReference`](@ref).
"""
function sph2xy(lon::Real, lat::Real, ref::EquirectangularReference; R::Real = 6371)
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

Compute spherical coordinates `[lon, lat]` [deg] from rectilinear coordinates `(x, y)` from with references `ref::EquirectangularReference`.

The units of `x` and `y` should be the same as the optional argument `R`, the radius of the Earth. The default is `R = 6371 km.`

Can be applied as `xy2sph(xy, ref; R = 6371)` where `xy` is a vector of vectors or a `Matrix`.

### Arguments

`x`: The x Cartesian coordinate in km (East/West).
`y`: The y Cartesian coordinate in km (North/South).
`ref`: An [`EquirectangularReference`](@ref).
"""
function xy2sph(x::Real, y::Real, ref::EquirectangularReference; R::Real = 6371)
    lon0, lat0 = (ref.lon0, ref.lat0)
    deg2rad = π/180
    rad2deg = 1/deg2rad

    lon = lon0 + rad2deg*x/(R*cos(lat0*deg2rad))
    lat = lat0 + rad2deg*y/R 

    return [lon, lat]
end

function xy2sph(xy::Vector{<:Vector{T}}, ref::EquirectangularReference; R::Real = 6371) where {T<:Real}
    lonlat = Matrix{T}(undef, length(xy), 2)
    
    for i = 1:length(xy)
        lonlat[i,:] = xy2sph(xy[i]...,ref, R = R)
    end

    return lonlat
end

function xy2sph(xy::Matrix{T}, ref::EquirectangularReference; R::Real = 6371) where {T<:Real}
    lonlat = Matrix{T}(undef, size(xy)...)
    
    for i = 1:size(xy, 1)
        lonlat[i,:] = xy2sph(xy[i,:]...,ref, R = R)
    end

    return lonlat
end