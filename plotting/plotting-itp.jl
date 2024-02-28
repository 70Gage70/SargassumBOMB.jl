"""
    check_land(land_itp; t, tol, limits, n_points)

Construct a plot of the land locations from `land_itp`.

### Arguments

- `land_itp`: A `InterpolatedField` which gives the interpolated location of the land.

### Optional Arguments

- `limits`: The `(lon_min, lon_max, lat_min, lat_max)` boundaries of the plot. Default `limits = (-180, 180, -90, 90)`.
- `n_points`: The number of points to use in each dimension of the plot, more gives higher resolution. Default `n_points = 1000`.
"""
function check_land(
    land_itp::InterpolatedField;
    limits::NTuple{4, Real} = (-180, 180, -90, 90),
    n_points::Integer = 1000)

    xs = range(start=limits[1], stop=limits[2], length=n_points)
    ys = range(start=limits[3], stop=limits[4], length=n_points)
    zs = [land_itp.fields[:land](sph2xy(x, y, land_itp.ref)...) for x in xs, y in ys]

    fig = default_fig()
    ax = geo_axis(fig[1, 1], title="Land", limits=limits)
    heatmap!(ax, xs, ys, zs)

    return fig
end


"""
    check_windwater(windwater_itp; t, tol, limits, n_points, u_name, v_name)

Construct a plot of the `u` and `v` components of `windwater_itp`.

### Arguments

- `windwater_itp`: A `Vector{InterpolatedField}` whose entries are wind/water vector fields. Each vector field is plotted in its 
    own row.
- `names`: A `Vector` of `String`s such that `names[i]` is the name of `windwater_itp[i]`.
- `time`: The time at which to plot the fields.

### Optional Arguments

- `limits`: The `(lon_min, lon_max, lat_min, lat_max)` boundaries of the plot. Default `limits = (-100, -50, 5, 35)`.
- `n_points`: The number of points to use in each dimension of the plot, more gives higher resolution. Default `n_points = 1000`.
- `u_name`: A `Symbol` giving the name of the field corresponding to the x component of the velocity. Default `u_name = :u`.
- `v_name`: A `Symbol` giving the name of the field corresponding to the x component of the velocity. Default `v_name = :v`.
"""
function check_windwater(
    windwater_itp::Vector{<:InterpolatedField},
    names::Vector{<:String},
    time::Real = 0.0;
    limits::NTuple{4, Real} = (-100, -50, 5, 35),
    n_points::Integer = 1000,
    u_name::Symbol = :u,
    v_name::Symbol = :v)

    @assert length(names) == length(windwater_itp)

    fig = default_fig()

    xs = range(start=limits[1], stop=limits[2], length = n_points)
    ys = range(start=limits[3], stop=limits[4], length = n_points)

    for i = 1:length(windwater_itp)
        itp = windwater_itp[i]
        us = [itp.fields[u_name](sph2xy(x, y)..., time) for x in xs, y in ys]
        vs = [itp.fields[v_name](sph2xy(x, y)..., time) for x in xs, y in ys]

        ax_u = geo_axis(fig[i, 1], title = "$(names[i]) u", limits = limits)
        heatmap!(ax_u, xs, ys, us)

        ax_v = geo_axis(fig[i, 2], title = "$(names[i]) v", limits = limits)
        heatmap!(ax_v, xs, ys, vs)
    end    

    return fig
end

"""
    vector_field_t!(axis, vf, t; fieldnames)
"""
function vector_field_t!(
    axis::Axis,
    vf::InterpolatedField,
    t::Real;
    fieldnames::Tuple{Symbol, Symbol} = (:u, :v))

    lims = axis.limits.val
    f(x, y) = Point2f(
        vf.fields[fieldnames[1]](sph2xy(x, y, vf.ref)..., t), 
        vf.fields[fieldnames[2]](sph2xy(x, y, vf.ref)..., t))
    
    streamplot!(axis, f, lims[1]..lims[2], lims[3]..lims[4])
end

"""
    scalar_field_t!(axis, sf, t; fieldnames, n_points)
"""
function scalar_field_t!(
    axis::Axis,
    sf::InterpolatedField,
    t::Real;
    fieldname::Symbol = :u,
    n_points::Integer = 100)

    lims = axis.limits.val

    xs = range(start = lims[1], stop = lims[2], length = n_points)
    ys = range(start = lims[3], stop = lims[4], length = n_points)
    zs = [sf.fields[fieldname](sph2xy(x, y, sf.ref)..., t) for x in xs, y in ys]
    
    heatmap!(axis, xs, ys, zs)
end