"""
    update_interpolant!(itp, itp_new)

Update (replace) `itp` with `itp_new`, where `itp_new` should be an [`InterpolatedField`](@ref) and \
`itp` should be one of 

- `WATER_ITP`
- `WIND_ITP`
- `STOKES_ITP`
- `WAVES_ITP`
- `NUTRIENTS_ITP`
- `TEMPERATURE_ITP`
- `LAND_ITP`
"""
function update_interpolant!(itp::Ref{InterpolatedField}, itp_new::InterpolatedField)
    itp.x = itp_new
    return nothing
end

"""
    limits(itp)

Return a the corners of the spatio-temporal box the interpolant is defined in.

If `itp` is time-dependent, return `(lon_min, lon_max, lat_min, lat_max, t_min, t_max)`.

If `itp` is time-independent, `(lon_min, lon_max, lat_min, lat_max)`.

### Examples

```julia-repl
(lon_min, lon_max, lat_min, lat_max, t_min, t_max) = limits(WATER_ITP)
```

```julia-repl
(lon_min, lon_max, lat_min, lat_max) = limits(LAND_ITP)
```
"""
function limits(itp)
	x = itp.x
	lon, lat = xy2sph(x.dims[:x], x.dims[:y]) .|> extrema
	lon = lon .|> x -> round(x, sigdigits = 4)
	lat = lat .|> x -> round(x, sigdigits = 4)

    if :t ∈ keys(x.dims)
	    tmin, tmax = extrema(x.dims[:t]) .|> time2datetime .|> x -> "$(monthname(x)), $(day(x)), $(year(x))"
	    return (lon[1], lon[2], lat[1], lat[2], tmin, tmax)
    else
        return (lon[1], lon[2], lat[1], lat[2])
    end
end

"""
    dim(itp, name)

Return the variable of `itp` whose value corresponds to the dimension indicated by `name`.

Use [`dims`](@ref) to see a list of possible values of `name`.

### Example

```julia-repl
x = dim(WATER_ITP, :x) # the x values defining interpolant knots
```
"""
dim(itp::Ref{InterpolatedField}, name::Symbol) = itp.x.dims[name] 

"""
    dims(itp)

Return the list of variable name/unit pairs of `itp`.

### Example

```julia-repl
dims(WATER_ITP)
```
"""
dims(itp::Ref{InterpolatedField}) = itp.x.dims_names

"""
    field(itp, name)

Return the sub-interpolant of `itp` whose value corresponds to the field indicated by `name`.

Use [`fields`](@ref) to see a list of possible values of `name`.

### Example

```julia-repl
v_x = field(WATER_ITP, :u) # the x component of the water velocity
v_x(1, 2, 3) # evaluate it at `(x, y, t) = (1, 2, 3)`.
```
"""
field(itp::Ref{InterpolatedField}, name::Symbol) = itp.x.fields[name] 

"""
    fields(itp)

Return the list of field name/unit pairs of `itp`.

### Example

```julia-repl
fields(WATER_ITP)
```
"""
fields(itp::Ref{InterpolatedField}) = itp.x.fields_names