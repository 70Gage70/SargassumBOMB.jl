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
    zs = [land_itp.fields[:land](sph2xy(x, y)...) for x in xs, y in ys]

    set_theme!(GEO_THEME())
    fig = Figure()
    ax = Axis(fig[1, 1], title="Land", limits=limits)
    heatmap!(ax, xs, ys, zs)

    return fig
end


"""
    check_itp(itp; title, time, limits, n_points, u_name, v_name)

Construct a plot of the `u_name` and `v_name` components of `itp`.

### Arguments

- `itp`: A [`InterpolatedField`](@ref) whose entries are wind/water vector fields.

### Optional Arguments

- `title`: A `String`s which labels the plots.
- `time`: The time at which to plot the fields. Default the minimum time.
- `limits`: The `(lon_min, lon_max, lat_min, lat_max)` boundaries of the plot. Default `limits = (-100, -50, 5, 35)`.
- `n_points`: The number of points to use in each dimension of the plot, more gives higher resolution. Default `n_points = 1000`.
- `u_name`: A `Symbol` giving the name of the field corresponding to the x component of the velocity. Default `u_name = :u`.
- `v_name`: A `Symbol` giving the name of the field corresponding to the y component of the velocity. Default `v_name = :v`.
"""
function check_itp(
    itp::InterpolatedField;
    title::String = "plot",
    time::Real = first(itp.dims[:t]),
    limits::NTuple{4, Real} = (-100, -50, 5, 35),
    n_points::Integer = 1000,
    u_name::Symbol = :u,
    v_name::Symbol = :v)

    set_theme!(GEO_THEME())
    fig = Figure()

    xs = range(start=limits[1], stop=limits[2], length = n_points)
    ys = range(start=limits[3], stop=limits[4], length = n_points)

    us = [itp.fields[u_name](sph2xy(x, y)..., time) for x in xs, y in ys]
    vs = [itp.fields[v_name](sph2xy(x, y)..., time) for x in xs, y in ys]

    ax_u = Axis(fig[1, 1], title = "$(title) u", limits = limits)
    heatmap!(ax_u, xs, ys, us)

    ax_v = Axis(fig[1, 2], title = "$(title) v", limits = limits)
    heatmap!(ax_v, xs, ys, vs)

    return fig
end

# """
#     vector_field_t!(axis, vf, t; fieldnames)
# """
# function vector_field_t!(
#     axis::Axis,
#     vf::InterpolatedField,
#     t::Real;
#     fieldnames::Tuple{Symbol, Symbol} = (:u, :v))

#     lims = axis.limits.val
#     f(x, y) = Point2f(
#         vf.fields[fieldnames[1]](sph2xy(x, y)..., t), 
#         vf.fields[fieldnames[2]](sph2xy(x, y)..., t))
    
#     streamplot!(axis, f, lims[1]..lims[2], lims[3]..lims[4])
# end

# """
#     scalar_field_t!(axis, sf, t; fieldnames, n_points)
# """
# function scalar_field_t!(
#     axis::Axis,
#     vf::InterpolatedField,
#     t::Real;
#     fieldname::Symbol = :u,
#     n_points::Integer = 100)

#     lims = axis.limits.val

#     xs = range(start = lims[1], stop = lims[2], length = n_points)
#     ys = range(start = lims[3], stop = lims[4], length = n_points)
#     zs = [vf.fields[fieldname](sph2xy(x, y)..., t) for x in xs, y in ys]
    
#     heatmap!(axis, xs, ys, zs)
# end