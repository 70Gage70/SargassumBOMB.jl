"""
    EquirectangularReference{T}

A container for the reference longitude and latitude of an equirectangular projection.

### Fields

- `lon0`: The standard longitude degrees (East/West).
- `lat0`: The standard latitude degrees (North/South).
- `R`: The radius of the Earth.
"""
struct EquirectangularReference{T<:Real}
    lon0::T
    lat0::T
    R::T
end

"""
    EquirectangularReference(; lon0, lat0)

Construct an `EquirectangularReference` with longitude `lon0` in degrees (East/West), latitude `lat0` in degrees (North/South) and with the radius of the Earth `R` (default 6371 km.)
"""
function EquirectangularReference(; lon0::Real, lat0::Real, R::Real = 6371.0)
    @assert -180.0 <= lon0 <= 180.0 "The longitude must be between -180 degrees and 180 degrees."
    @assert -90 <= lat0 <= 90 "The latitude must be between -90 degrees and 90 degrees."
    @assert R > 0 "The radius of the Earth should be positive."

    lon0, lat0, R = promote(float(lon0), float(lat0), float(R))

    return EquirectangularReference(lon0, lat0, R)
end

"""
    sph2xy(lon, lat, ref)

Compute planar coordinates `[x, y]` from spherical coordinates `(lon, lat)` [deg] with reference `ref::EquirectangularReference`.

The units of `x` and `y` the same as `ref.R`.

Can be applied as `sph2xy(lon_range, lat_range, ref)` where `lon_range` and `lat_range` are `AbstractRange`. Returns `(x_range, y_range)`.

### Arguments

- `lon`: Longitude in degrees (East/West).
- `lat`: Latitude in degrees (North/South).
- `ref`: An [`EquirectangularReference`](@ref).
"""
function sph2xy(lon::Real, lat::Real, ref::EquirectangularReference)
    @assert -180.0 <= lon <= 180.0 "The longitude must be between -180 degrees and 180 degrees."
    @assert -90 <= lat <= 90 "The latitude must be between -90 degrees and 90 degrees."

    lon0, lat0, R = (ref.lon0, ref.lat0, ref.R)
    deg2rad = π/180

    x = R*(lon - lon0)*deg2rad*cos(lat0*deg2rad)
    y = R*(lat - lat0)*deg2rad

    return [x, y]
end

function sph2xy(lon_range::AbstractRange, lat_range::AbstractRange, ref::EquirectangularReference)
    # uses the fact that the translation between eqr and spherical is linear
    lonmin, latmin = sph2xy(first(lon_range), first(lat_range), ref)
    lonmax, latmax = sph2xy(last(lon_range), last(lat_range), ref)

    return  (
            range(start = lonmin, length = length(lon_range), stop = lonmax), 
            range(start = latmin, length = length(lat_range), stop = latmax)
            )
end

"""
    sph2xy(gridded_field, ref; lon_name = :lon, lat_name = :lat, x_name = :x, y_name = :y)

Transform the `lon`, `lat` variables in `gridded_field` to `x`, `y` variables with reference `ref::EquirectangularReference` and
return a new `GriddedField`.

The units of `x` and `y` are the same as `ref.R`.

The `lon` and `lat` variables in `gridded_field` should have names `lat_name` (default `:lat`) and `lon_name` (default `:lon`).

The returned `GriddedField` has its `ref` field updated to the input argument `ref` and its variable names updated to `x_name` (default `:x`)
and `y_name` (default `:y`).
"""
function sph2xy(
    gridded_field::GriddedField, 
    ref::EquirectangularReference;
    lon_name::Symbol = :lon, 
    lat_name::Symbol = :lat, 
    x_name::Symbol = :x, 
    y_name::Symbol = :y)

    lon = gridded_field.vars[lon_name]
    lat = gridded_field.vars[lat_name]
    x, y = sph2xy(lon, lat, ref)

    new_var_names = Vector{Symbol}()
    new_vars = typeof(lon)()

    for var_name in gridded_field.var_names
        if var_name == lon_name
            new_vars[x_name] = x
            push!(new_var_names, x_name)
        elseif var_name == lat_name
            new_vars[y_name] = y
            push!(new_var_names, y_name)
        else
            new_vars[var_name] = gridded_field.vars[var_name]
            push!(new_var_names, var_name)
        end
    end

    return GriddedField(new_var_names, new_vars, gridded_field.fields, gridded_field.vars_units, gridded_field.fields_units, gridded_field.time_start, ref)
end

"""
    xy2sph(x, y, ref)

Compute spherical coordinates `[lon, lat]` [deg] from rectilinear coordinates `(x, y)` from with references `ref::EquirectangularReference`.

The units of `x` and `y` should be the same as `ref.R`.

Can be applied as `xy2sph(xy, ref)` where `xy` is a `Vector{Vector{T}}` or a `Matrix`. Returns a `Matrix`.

Can be applied as `xy2sph(x_range, y_range, ref)` where `x_range` and `y_range` are `AbstractRange`. Returns `(lon_range, lat_range)`.

### Arguments

- `x`: The x Cartesian coordinate in km (East/West).
- `y`: The y Cartesian coordinate in km (North/South).
- `ref`: An [`EquirectangularReference`](@ref).
"""
function xy2sph(x::Real, y::Real, ref::EquirectangularReference)
    lon0, lat0, R = (ref.lon0, ref.lat0, ref.R)
    deg2rad = π/180
    rad2deg = 1/deg2rad

    lon = lon0 + rad2deg*x/(R*cos(lat0*deg2rad))
    lat = lat0 + rad2deg*y/R 

    return [lon, lat]
end

function xy2sph(xy::Vector{<:Vector{T}}, ref::EquirectangularReference) where {T<:Real}
    lonlat = zeros(length(xy), 2) 
    
    for i = 1:length(xy)
        lonlat[i,:] = xy2sph(xy[i]..., ref)
    end

    return lonlat
end

function xy2sph(xy::Matrix{T}, ref::EquirectangularReference) where {T<:Real}
    lonlat = zeros(size(xy)...)
    
    for i = 1:size(xy, 1)
        lonlat[i,:] = xy2sph(xy[i,:]..., ref)
    end

    return lonlat
end

function xy2sph(x_range::AbstractRange, y_range::AbstractRange, ref::EquirectangularReference)
    # uses the fact that the translation between eqr and spherical is linear
    xmin, ymin = xy2sph(first(x_range), first(y_range), ref)
    xmax, ymax = xy2sph(last(x_range), last(y_range), ref)

    return  (
            range(start = xmin, length = length(x_range), stop = xmax), 
            range(start = ymin, length = length(y_range), stop = ymax)
            )
end