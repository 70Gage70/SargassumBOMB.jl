using MAT
using Interpolations

include("itp-helpers.jl")
include(joinpath(@__DIR__, "..", "src", "coordinates.jl"))

############################################################################################

###########################
########################### GRIDS
###########################

"""
    AbstractVectorFieldGrid

The abstract type for gridded vector fields.
"""
abstract type AbstractVectorFieldGrid end

"""
    mutable struct VectorField2DGridSPH{T}

A container for a two-dimensional vector field defined at each point on a linearly spaced grid of longitude, latitude and time.

### Fields

- `lon`: The range of longitudes (East/West) over which the grid is defined.
- `lat`: The range of latitudes (North/South) over which the grid is defined.
- `time0`: The initial `DateTime` at which the vector field is defined.
- `time`: The range of times since `time0` over which the grid is defined.
- `u`: The x component of the vector field, conventionally u[lon, lat, time].
- `v`: The y component of the vector field, conventionally v[lon, lat, time].
"""
mutable struct VectorField2DGridSPH{A<:AbstractRange, T<:Real} <: AbstractVectorFieldGrid
    lon::A
    lat::A
    time0::DateTime
    time::A
    u::Array{T,3}
    v::Array{T,3}
end

"""
    VectorField2DGridSPH(infile; kwargs...)

Create a `VectorField2DGridSPH` object from the file `infile`. 

`infile` must be of the form `"filename.mat"`.

The longitudes, latitudes and times must all be defined on linearly spaced grids. In particular, the times should be given as
days in Rata Die format. See [`rata2datetime_minute`](@ref).

### Optional Arguments

- `lon_alias`: The name of the variable in `infile` which stores the longitudes. Default `"lon"`.
- `lat_alias`: The name of the variable in `infile` which stores the latitudes. Default `"lat"`.
- `time_alias`: The name of the variable in `infile` which stores the times. Default `"t"`.
- `u_alias`: The name of the variable in `infile` which stores the x component of the vector field. Default `"u"`.
- `v_alias`: The name of the variable in `infile` which stores the y component of the vector field. Default `"v"`.
- `NaN_replacement`: Any `NaN`s in `u` and `v` will be replaced by this. Default `0.0`.
- `lon_lat_time_order`: A permutation that will be applied to `u` and `v`. Use if `u`, `v` do not follow the convention [lon, lat, time]. Default `[1, 2, 3]`.
"""
function VectorField2DGridSPH(
    infile::String;
    lon_alias::String = "lon",
    lat_alias::String = "lat",
    time_alias::String = "t",
    u_alias::String = "u",
    v_alias::String = "v",
    NaN_replacement::Any = 0.0,
    lon_lat_time_order::Vector{<:Integer} = [1, 2, 3])

    ### loading data

    # ensure that the provided file is a mat file
    extension = infile[findlast(==('.'), infile)+1:end]
    @assert extension in ["mat"] "Require a .mat file."

    data = matopen(infile)
    data_keys = collect(keys(data))

    # ensure that all the relevant variables are contained in the given file
    @assert lon_alias in data_keys "$(lon_alias) not in $(infile)"
    @assert lat_alias in data_keys "$(lat_alias) not in $(infile)"
    @assert time_alias in data_keys "$(time_alias) not in $(infile)"
    @assert u_alias in data_keys "$(u_alias) not in $(infile)"
    @assert v_alias in data_keys "$(v_alias) not in $(infile)"

    # lon, lat, time should all be vectors, as opposed to N x 1 matrices.
    lon, lat, time = map(vec, read(data, lon_alias, lat_alias, time_alias));
    u, v = read(data, u_alias, v_alias);
    close(data)

    ### checking dimension consistency
    
    u = permutedims(u, lon_lat_time_order)
    v = permutedims(v, lon_lat_time_order)  
    
    u_size = size(u)
    v_size = size(v)
    
    @assert u_size == v_size "The x and y components of the vector field must have the same dimensions. Got size(u) = $(u_size) and size(v) = $(v_size)."
    
    # ensure that u and v follow the convention [lon, lat, time]
    @assert length(lon) == u_size[1] "The number of entries in the longitude [$(length(lon))] must match the first dimension of the vector components [$(u_size[1])]. Use the kwarg `lon_lat_time_order` to supply a permutation if the vector components are out of order."

    @assert length(lat) == u_size[2] "The number of entries in the latitude [$(length(lat))] must match the second dimension of the vector components [$(u_size[2])]. Use the kwarg `lon_lat_time_order` to supply a permutation if the vector components are out of order."

    @assert length(time) == u_size[3] "The number of entries in the time [$(length(time))] must match the third dimension of the vector components [$(u_size[3])]. Use the kwarg `lon_lat_time_order` to supply a permutation if the vector components are out of order."

    ### cleaning data

    # convert linearly spaced vector to a range
    lon = reduce_vector_to_range(lon)
    lat = reduce_vector_to_range(lat)
    time = reduce_vector_to_range(time)

    # ensure that the ranges are increasing (if they are decreasing, simply reverse the appropriate dimension of u and v and corresponding variable)
    if step(lon) < 0
        reverse!(u, dims = 1)
        reverse!(v, dims = 1)
        lon = reverse(lon)
    end

    if step(lat) < 0
        reverse!(u, dims = 2)
        reverse!(v, dims = 2)
        lat = reverse(lat)
    end

    if step(time) < 0
        reverse!(u, dims = 3)
        reverse!(v, dims = 3)
        time = reverse(time)
    end

    # ensure that longitudes and latitudes are in the allowed earthly range
    @assert (-180.0 <= first(lon) <= 180.0) && (-180.0 <= last(lon) <= 180.0) "The given longitudes are not between -180 degrees and 180 degrees."
    @assert (-90.0 <= first(lat) <= 90.0) && (-90.0 <= last(lat) <= 90.0) "The given latitudes are not between -90 degrees and 90 degrees."

    # NaN replacement
    u[isnan.(u)] .= NaN_replacement
    v[isnan.(v)] .= NaN_replacement

    # rescale time
    time0 = first(time)
    time = time .- time0
    time0 = rata2datetime_minute(time0)

    # ensure data are of matching types
    lon, lat, time = promote(lon, lat, time)
    u, v = promote(u, v)

    return VectorField2DGridSPH(lon, lat, time0, time, u, v)
end

"""
    AbstractScalarFieldGrid

The abstract type for gridded scalar fields.
"""
abstract type AbstractScalarFieldGrid end

"""
    mutable struct ScalarField2DGridSPH{T}

A container for a two-dimensional scalar field defined at each point on a linearly spaced grid of longitude, latitude and time.

### Fields

- `lon`: The range of longitudes (East/West) over which the grid is defined.
- `lat`: The range of latitudes (North/South) over which the grid is defined.
- `time0`: The initial `DateTime` at which the vector field is defined.
- `time`: The range of times since `time0` over which the grid is defined.
- `u`: The value of the scalar field, conventionally u[lon, lat, time].
"""
mutable struct ScalarField2DGridSPH{A<:AbstractRange, T<:Real} <: AbstractScalarFieldGrid
    lon::A
    lat::A
    time0::DateTime
    time::A
    u::Array{T,3}
end

"""
    ScalarField2DGridSPH(infile; kwargs...)

Create a `ScalarField2DGridSPH` object from the file `infile`. 

`infile` must be of the form `"filename.mat"`.

The longitudes, latitudes and times must all be defined on linearly spaced grids. In particular, the times should be given as
days in Rata Die format. See [`rata2datetime_minute`](@ref).

### Optional Arguments

- `lon_alias`: The name of the variable in `infile` which stores the longitudes. Default `"lon"`.
- `lat_alias`: The name of the variable in `infile` which stores the latitudes. Default `"lat"`.
- `time_alias`: The name of the variable in `infile` which stores the times. Default `"t"`.
- `u_alias`: The name of the variable in `infile` which stores the x component of the vector field. Default `"u"`.
- `NaN_replacement`: Any `NaN`s in `u` and `v` will be replaced by this. Default `0.0`.
- `lon_lat_time_order`: A permutation that will be applied to `u`. Use if `u` does not follow the convention [lon, lat, time]. Default `[1, 2, 3]`.
"""
function ScalarField2DGridSPH(
    infile::String;
    lon_alias::String = "lon",
    lat_alias::String = "lat",
    time_alias::String = "t",
    u_alias::String = "u",
    NaN_replacement::Any = 0.0,
    lon_lat_time_order::Vector{<:Integer} = [1, 2, 3])

    vf = VectorField2DGridSPH(
        infile,
        lon_alias = lon_alias,
        lat_alias = lat_alias,
        time_alias = time_alias,
        u_alias = u_alias,
        v_alias = u_alias, # just set v = u
        NaN_replacement = NaN_replacement,
        lon_lat_time_order = lon_lat_time_order
    )

    return ScalarField2DGridSPH(vf.lon, vf.lat, vf.time0, vf.time, vf.u)
end

###########################
########################### INTERPOLANTS
###########################

###########################
########################### VECTOR FIELDS
###########################

"""
    const DEFAULT_INTERPOLATION

The default interpolation method, `BSpline(Cubic(Line(OnGrid())))`.
"""
const DEFAULT_INTERPOLATION = BSpline(Cubic(Interpolations.Line(OnGrid())))

"""
    interpolate_field(x, y, t, u; interpolant_type)

Construct an interpolant for the field `u[x, y, t]`. Outside of the range defined by `x`, `y` and `t` defaults to 0.0 via extrapolation.

### Optional Arguments
- `interpolant_type`: The type of the interpolant used at the interpolation stage. Default to [DEFAULT_INTERPOLATION](@ref).
"""
function interpolate_field(
    x::AbstractRange, 
    y::AbstractRange, 
    t::AbstractRange, 
    u::AbstractArray; 
    interpolant_type = DEFAULT_INTERPOLATION)

    return extrapolate(scale(Interpolations.interpolate(u, interpolant_type), x, y, t), 0.0)
end

"""
    AbstractVectorFieldInterpolant

The abstract type for interpolated vector fields.
"""
abstract type AbstractVectorFieldInterpolant end

"""
    struct VectorField2DInterpolantSPH{T, I}

A container for a two-dimensional vector field interpolated over spherical longitude, latitude and time coordinates.

### Fields

- `lon`: The range of longitudes (East/West) over which the grid is defined.
- `lat`: The range of latitudes (North/South) over which the grid is defined.
- `time0`: The initial `DateTime` at which the vector field is defined.
- `time`: The range of times since `time0` over which the grid is defined.
- `u`: An interpolant for the x component of the vector field. Call as `u(lon, lat, t)`.
- `v`: An interpolant for the y component of the vector field. Call as `v(lon, lat, t)`.
"""
struct VectorField2DInterpolantSPH{R<:AbstractRange, I<:AbstractInterpolation} <: AbstractVectorFieldInterpolant
    lon::R
    lat::R
    time0::DateTime
    time::R
    u::I
    v::I
end

"""
    VectorField2DInterpolantSPH(vf_grid; interpolant_type)

Construct an `VectorField2DInterpolantSPH` from a [`VectorField2DGridSPH`](@ref) using [`interpolate_field`](@ref).
"""
function VectorField2DInterpolantSPH(vf_grid::VectorField2DGridSPH)
    lon, lat, time, time0, u, v = (vf_grid.lon, vf_grid.lat, vf_grid.time, vf_grid.time0, vf_grid.u, vf_grid.v)

    u_itp = interpolate_field(lon, lat, time, u, interpolant_type = interpolant_type)
    v_itp = interpolate_field(lon, lat, time, v, interpolant_type = interpolant_type)

    return VectorField2DInterpolantSPH(lon, lat, time0, time, u_itp, v_itp)
end

"""
    struct VectorField2DInterpolantEQR{T, I}

A container for a two-dimensional vector field interpolated over equirectangular x, y and time coordinates.

### Fields

- `ref`: The `EquirectangularReference` with whcih the projection is defined.
- `x`: The range of x values over which the grid is defined.
- `y`: The range of y values over which the grid is defined.
- `time0`: The initial `DateTime` at which the vector field is defined.
- `time`: The range of times since `time0` over which the grid is defined.
- `u`: An interpolant for the x component of the vector field. Call as `u(x, y, t)`.
- `v`: An interpolant for the y component of the vector field. Call as `v(x, y, t)`.
"""
struct VectorField2DInterpolantEQR{T<:Real, R<:AbstractRange, I<:AbstractInterpolation} <: AbstractVectorFieldInterpolant
    ref::EquirectangularReference{T}
    x::R
    y::R
    time0::DateTime
    time::R
    u::I
    v::I
end

"""
    VectorField2DInterpolantEQR(vf_grid, ref; interpolant_type)

Construct an `VectorField2DInterpolantEQR` from a [`VectorField2DGridSPH`](@ref) and `ref::EquirectangularReference` using [`interpolate_field`](@ref).

### Optional Arguments

- `interpolant_type`: The type of interpolation to use from the `Interpolations` module. Default cubic B-spline.
"""
function VectorField2DInterpolantEQR(
    vf_grid::VectorField2DGridSPH,
    ref::EquirectangularReference; 
    interpolant_type = DEFAULT_INTERPOLATION
    )

    lon, lat, time, time0, u, v = (vf_grid.lon, vf_grid.lat, vf_grid.time, vf_grid.time0, vf_grid.u, vf_grid.v)

    x, y = sph2xy(lon, lat, ref)

    u_itp = interpolate_field(x, y, time, u, interpolant_type = interpolant_type)
    v_itp = interpolate_field(x, y, time, v, interpolant_type = interpolant_type)

    return VectorField2DInterpolantEQR(ref, x, y, time0, time, u_itp, v_itp)
end

###########################
########################### SCALAR FIELDS
###########################

"""
    AbstractScalarFieldInterpolant

The abstract type for interpolated scalar fields.
"""
abstract type AbstractScalarFieldInterpolant end

"""
    struct ScalarField2DInterpolantSPH{T, I}

A container for a two-dimensional scalar field interpolated over spherical longitude, latitude and time coordinates.

### Fields

- `lon`: The range of longitudes (East/West) over which the grid is defined.
- `lat`: The range of latitudes (North/South) over which the grid is defined.
- `time0`: The initial `DateTime` at which the scalar field is defined.
- `time`: The range of times since `time0` over which the grid is defined.
- `u`: An interpolant for the scalar field. Call as `u(lon, lat, t)`.
"""
struct ScalarField2DInterpolantSPH{R<:AbstractRange, I<:AbstractInterpolation} <: AbstractScalarFieldInterpolant
    lon::R
    lat::R
    time0::DateTime
    time::R
    u::I
end

"""
    ScalarField2DInterpolantSPH(sf_grid; interpolant_type)

Construct an `ScalarField2DInterpolantSPH` from a [`ScalarField2DGridSPH`](@ref) using [`interpolate_field`](@ref).
"""
function ScalarField2DInterpolantSPH(sf_grid::ScalarField2DGridSPH)
    lon, lat, time, time0, u = (sf_grid.lon, sf_grid.lat, sf_grid.time, sf_grid.time0, sf_grid.u)

    u_itp = interpolate_field(lon, lat, time, u, interpolant_type = interpolant_type)

    return ScalarField2DInterpolantSPH(lon, lat, time0, time, u_itp)
end

"""
    struct ScalarField2DInterpolantEQR{T, I}

A container for a two-dimensional scalar field interpolated over equirectangular x, y and time coordinates.

### Fields

- `ref`: The `EquirectangularReference` with whcih the projection is defined.
- `x`: The range of x values over which the grid is defined.
- `y`: The range of y values over which the grid is defined.
- `time0`: The initial `DateTime` at which the scalar field is defined.
- `time`: The range of times since `time0` over which the grid is defined.
- `u`: An interpolant for the scalar field. Call as `u(x, y, t)`.
"""
struct ScalarField2DInterpolantEQR{T<:Real, R<:AbstractRange, I<:AbstractInterpolation} <: AbstractScalarFieldInterpolant
    ref::EquirectangularReference{T}
    x::R
    y::R
    time0::DateTime
    time::R
    u::I
end

"""
    ScalarField2DInterpolantEQR(vf_grid, ref; interpolant_type)

Construct an `ScalarField2DInterpolantEQR` from a [`ScalarField2DGridSPH`](@ref) and `ref::EquirectangularReference` using [`interpolate_field`](@ref).

### Optional Arguments

- `interpolant_type`: The type of interpolation to use from the `Interpolations` module. Default cubic B-spline.
"""
function ScalarField2DInterpolantEQR(
    sf_grid::ScalarField2DGridSPH,
    ref::EquirectangularReference; 
    interpolant_type = DEFAULT_INTERPOLATION
    )

    lon, lat, time, time0, u = (sf_grid.lon, sf_grid.lat, sf_grid.time, sf_grid.time0, sf_grid.u)

    x, y = sph2xy(lon, lat, ref)

    u_itp = interpolate_field(x, y, time, u, interpolant_type = interpolant_type)

    return ScalarField2DInterpolantEQR(ref, x, y, time0, time, u_itp)
end

###########################
########################### STATIC FIELDS
###########################

"""
    AbstractStaticFieldInterpolant

The abstract type for interpolated static fields.
"""
abstract type AbstractStaticFieldInterpolant end

"""
    struct StaticField2DInterpolantEQR{T, I}

A container for a two-dimensional static field interpolated over equirectangular x, y coordinates.

### Fields

- `ref`: The `EquirectangularReference` with whcih the projection is defined.
- `x`: The range of x values over which the grid is defined.
- `y`: The range of y values over which the grid is defined.
- `u`: An interpolant for the Static field. Call as `u(x, y)`.
"""
struct StaticField2DInterpolantEQR{T<:Real, R<:AbstractRange, I<:AbstractInterpolation} <: AbstractStaticFieldInterpolant
    ref::EquirectangularReference{T}
    x::R
    y::R
    u::I
end

"""
    StaticField2DInterpolantEQR(vf_grid, ref; interpolant_type)

Construct an `StaticField2DInterpolantEQR`.

### Optional Arguments

- `interpolant_type`: The type of interpolation to use from the `Interpolations` module. Default nearest-neighbor.
"""
function StaticField2DInterpolantEQR(
    lon::AbstractRange,
    lat::AbstractRange,
    data::AbstractArray,
    ref::EquirectangularReference; 
    interpolant_type = BSpline(Constant())
    )

    x, y = sph2xy(lon, lat, ref)

    u_itp = extrapolate(scale(interpolate(data, interpolant_type), x, y), 0)

    return StaticField2DInterpolantEQR(ref, x, y, u_itp)
end