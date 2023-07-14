using MAT
using Interpolations
using LinearAlgebra: ⋅
import LinearAlgebra.BLAS.nrm2 as norm

include(joinpath(@__DIR__, "..", "helpers.jl"))
include(joinpath(@__DIR__, "..", "parameters.jl"))

##############################################

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
mutable struct VectorField2DGridSPH{T<:Real} <: AbstractVectorFieldGrid
    lon::AbstractRange{T}
    lat::AbstractRange{T}
    time0::DateTime
    time::AbstractRange{T}
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
    interpolate_field(x, y, t, u; interpolant_type)

Construct an interpolant for the field `u[x, y, t]`. Outside of the range defined by `x`, `y` and `t` defaults to 0.0 via extrapolation.

### Optional Arguments
- `interpolant_type`: The type of the interpolant used at the interpolation stage. Default a cubic spline, `BSpline(Cubic(Line(OnGrid()))))`.
"""
function interpolate_field(
    x::AbstractRange, 
    y::AbstractRange, 
    t::AbstractRange, 
    u::AbstractArray; 
    interpolant_type = BSpline(Cubic(Line(OnGrid()))))

    return extrapolate(scale(interpolate(u, interpolant_type), x, y, t), 0.0)
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
struct VectorField2DInterpolantSPH{T<:Real, I<:Interpolations.InterpolationType} <: AbstractVectorFieldInterpolant
    lon::AbstractRange{T}
    lat::AbstractRange{T}
    time0::DateTime
    time::AbstractRange{T}
    u::AbstractInterpolation{T, 3, I}
    v::AbstractInterpolation{T, 3, I}
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
struct VectorField2DInterpolantEQR{T<:Real, I<:Interpolations.InterpolationType} <: AbstractVectorFieldInterpolant
    ref::EquirectangularReference{T}
    x::AbstractRange{T}
    y::AbstractRange{T}
    time0::DateTime
    time::AbstractRange{T}
    u::AbstractInterpolation{T, 3, I}
    v::AbstractInterpolation{T, 3, I}
end

"""
    VectorField2DInterpolantEQR(vf_grid, ref; interpolant_type, R)

Construct an `VectorField2DInterpolantEQR` from a [`VectorField2DGridSPH`](@ref) and `ref::EquirectangularReference` using [`interpolate_field`](@ref).

### Optional Arguments

- `interpolant_type`: The type of interpolation to use from the `Interpolations` module. Default cubic B-spline.
- `R`: The radius of the Earth. The units of the `x` and `y` coordinates are the same as the units of `R`. Default 6371 km.
"""
function VectorField2DInterpolantEQR(
    vf_grid::VectorField2DGridSPH,
    ref::EquirectangularReference; 
    interpolant_type = BSpline(Cubic(Line(OnGrid()))),
    R = 6371
    )

    lon, lat, time, time0, u, v = (vf_grid.lon, vf_grid.lat, vf_grid.time, vf_grid.time0, vf_grid.u, vf_grid.v)

    x, y = sph2xy(lon, lat, ref, R = R)

    u_itp = interpolate_field(x, y, time, u, interpolant_type = interpolant_type)
    v_itp = interpolate_field(x, y, time, v, interpolant_type = interpolant_type)

    return VectorField2DInterpolantEQR(ref, x, y, time0, time, u_itp, v_itp)
end

"""
    MaterialDerivativeX(vector_field, x, y, t)

Compute the material derivative of the x component of the vector field `vector_field` at coordinates `(x, y, t)`.
"""
function MaterialDerivativeX(vector_field::VectorField2DInterpolantEQR, x::Real, y::Real, t::Real)
    return gradient(vector_field.u, x, y, t) ⋅ [vector_field.u(x, y, t), vector_field.v(x, y, t), 1.0]
end

"""
    MaterialDerivativeY(vector_field, x, y, t)

Compute the material derivative of the y component of the vector field `vector_field` at coordinates `(x, y, t)`.
"""
function MaterialDerivativeY(vector_field::VectorField2DInterpolantEQR, x::Real, y::Real, t::Real)
    return gradient(vector_field.v, x, y, t) ⋅ [vector_field.u(x, y, t), vector_field.v(x, y, t), 1.0]
end

"""
    Vorticity(vector_field, x, y, t)

Compute the vorticity of the vector field `vector_field` at coordinates `(x, y, t)`.
"""
function Vorticity(vector_field::VectorField2DInterpolantEQR, x::Real, y::Real, t::Real)
    return gradient(vector_field.v, x, y, t)[1] - gradient(vector_field.u, x, y, t)[2]
end