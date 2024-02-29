const UFUL = Unitful.Unitlike

"""
    const UNITS

A dictionary mapping dimension names to the `Unitful.Unitlike` that measures it.
"""
const UNITS = Dict{String, UFUL}(
    "distance" => u"km",
    "speed" => u"km/d",
    "wave_height" => u"m",
    "temperature" => u"°C",
    "concentration" => u"mmol/m^3",
    "degrees" => u"°",
    "time" => u"d",
    "none" => NoUnits
)

"""
    const EARTH_RADIUS

The radius of the Earth, equal to 6371 km.
"""
const EARTH_RADIUS = 6371.0 * u"km"

"""
    EquirectangularReference{T, U}

A container for the reference longitude and latitude of an equirectangular projection.

### Fields

- `lon0`: The standard longitude degrees (East/West).
- `lat0`: The standard latitude degrees (North/South).
- `R`: The radius of the Earth. The units of this quantity are the units of the equirectangular coordinates.

### Constructor

Use `EquirectangularReference(; lon0 = -75.0, lat0 = 10.0, units = UNITS["distance"])`.
"""
struct EquirectangularReference{T<:AbstractFloat, U<:Unitful.AbstractQuantity}
    lon0::T
    lat0::T
    R::U

    function EquirectangularReference(; lon0::Real = -75.0, lat0::Real = 10.0, units::Unitful.Unitlike = UNITS["distance"])
        @assert -180.0 <= lon0 <= 180.0 "The longitude must be between -180 degrees and 180 degrees."
        @assert -90 <= lat0 <= 90 "The latitude must be between -90 degrees and 90 degrees."
    
        lon0, lat0 = promote(float(lon0), float(lat0))
        R = uconvert(units, EARTH_RADIUS)
    
        return new{typeof(lon0), typeof(R)}(lon0, lat0, R)
    end
end

function Base.show(io::IO, x::EquirectangularReference)
    print(io, "EquirectangularReference[lon0 = $(x.lon0), lat0 = $(x.lat0), R = $(x.R)]")
end

"""
    const EQR
    
The [`EquirectangularReference`](@ref) used during all conversions. 

### Defaults

- `lon0`: -75.0 degrees
- `lat0`: 10.0 degrees
- `R`: 6731 km 
"""
const EQR = Ref{EquirectangularReference}(EquirectangularReference())

"""
    sph2xy(lon, lat)

Compute planar coordinates `[x, y]` from spherical coordinates `(lon, lat)` [deg].

The units of `x` and `y` the same as `EQR.R`.

Can be applied as `sph2xy(lon_range, lat_range)` where `lon_range` and `lat_range` are `AbstractRange`. Returns `(x_range, y_range)`.

Can be applied as `sph2xy(lon_lat)` where `lon_lat` is an `N x 2` `Matrix` or a `Vector` of length `2N` with 
entries of the form `[lon1, lat1, lon2, lat2 ... lonN, latN]`. Returns a result in the same shape as the input.

### Arguments

- `lon`: Longitude in degrees (East/West).
- `lat`: Latitude in degrees (North/South).
"""
function sph2xy(lon::Real, lat::Real)
    @assert -180.0 <= lon <= 180.0 "The longitude must be between -180 degrees and 180 degrees."
    @assert -90 <= lat <= 90 "The latitude must be between -90 degrees and 90 degrees."

    lon0, lat0, R = (EQR.x.lon0, EQR.x.lat0, EQR.x.R.val)
    deg2rad = π/180

    x = R*(lon - lon0)*deg2rad*cos(lat0*deg2rad)
    y = R*(lat - lat0)*deg2rad

    return [x, y]
end

function sph2xy(lon_range::AbstractRange, lat_range::AbstractRange)
    # uses the fact that the translation between eqr and spherical is linear
    lonmin, latmin = sph2xy(first(lon_range), first(lat_range))
    lonmax, latmax = sph2xy(last(lon_range), last(lat_range))

    return  (
            range(start = lonmin, length = length(lon_range), stop = lonmax), 
            range(start = latmin, length = length(lat_range), stop = latmax)
            )
end

function sph2xy(lon_lat::Matrix{T}) where {T<:Real}
    @assert size(lon_lat, 2) == 2 "lon_lat should be an `N x 2` matrix"
    xy = zeros(T, size(lon_lat))
    
    for i = 1:size(lon_lat, 1)
        xy[i,:] .= sph2xy(lon_lat[i,1], lon_lat[i,2])
    end

    return xy
end

function sph2xy(lon_lat::Vector{T}) where {T<:Real}
    @assert iseven(length(lon_lat)) "lon_lat should be of the form `[lon1, lat1, lon2, lat2 ... lon3, lat3]`."

    xy = zeros(T, length(lon_lat))
    
    for i = 1:2:length(lon_lat)
        xy[i:i+1] .= sph2xy(lon_lat[i], lon_lat[i + 1])
    end

    return xy
end

"""
    xy2sph(x, y)

Compute spherical coordinates `[lon, lat]` [deg] from rectilinear coordinates `(x, y)`.

The units of `x` and `y` should be the same as `EQR.R`.

Can be applied as `xy2sph(xy)` where `xy` is a `Vector{Vector{T}}` or an `N x 2` `Matrix`. Returns an `N x 2` `Matrix` in these cases.

Can be applied as `xy2sph(xy)` where `xy` a `Vector` of length `2N` with 
entries of the form `[lon1, lat1, lon2, lat2 ... lonN, latN]`. Returns a result in the same shape as the input.

Can be applied as `xy2sph(x_range, y_range)` where `x_range` and `y_range` are `AbstractRange`. Returns `(lon_range, lat_range)`.

### Arguments

- `x`: The x Cartesian coordinate in km (East/West).
- `y`: The y Cartesian coordinate in km (North/South).
"""
function xy2sph(x::Real, y::Real)
    lon0, lat0, R = (EQR.x.lon0, EQR.x.lat0, EQR.x.R.val)
    deg2rad = π/180
    rad2deg = 1/deg2rad

    lon = lon0 + rad2deg*x/(R*cos(lat0*deg2rad))
    lat = lat0 + rad2deg*y/R 

    return [lon, lat]
end

function xy2sph(xy::Vector{<:Vector{T}}) where {T<:Real}
    lonlat = zeros(T, length(xy), 2) 
    
    for i = 1:length(xy)
        lonlat[i,:] = xy2sph(xy[i][1], xy[i][2])
    end

    return lonlat
end

function xy2sph(xy::Matrix{T}) where {T<:Real}
    @assert size(xy, 2) == 2 "xy should be an `N x 2` matrix"
    lonlat = zeros(T, size(xy))
    
    for i = 1:size(xy, 1)
        lonlat[i,:] = xy2sph(xy[i,1], xy[i,2])
    end

    return lonlat
end

function xy2sph(x_range::AbstractRange, y_range::AbstractRange)
    # uses the fact that the translation between eqr and spherical is linear
    xmin, ymin = xy2sph(first(x_range), first(y_range))
    xmax, ymax = xy2sph(last(x_range), last(y_range))

    return  (
            range(start = xmin, length = length(x_range), stop = xmax), 
            range(start = ymin, length = length(y_range), stop = ymax)
            )
end

function xy2sph(xy::Vector{T}) where {T<:Real}
    @assert iseven(length(xy)) "xy should be of the form `[x1, y1, x2, y2 ... x3, y3]`."

    lon_lat = zeros(T, length(xy))
    
    for i = 1:2:length(xy)
        lon_lat[i:i+1] .= xy2sph(xy[i], xy[i + 1])
    end

    return lon_lat
end