include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/statistic-methods.jl"))

################################################################################

"""
    check_windwater(windwater_itp; t, tol, limits, n_points, u_name, v_name)

Construct a plot of the `u` and `v` components of `windwater_itp`.

### Arguments

- `windwater_itp`: A `InterpolatedField` for either wind or water speeds.
- `time`: The time at which to plot the fields.

### Optional Arguments

- `limits`: The `(lon_min, lon_max, lat_min, lat_max)` boundaries of the plot. Default `limits = (-100, -50, 5, 35)`.
- `n_points`: The number of points to use in each dimension of the plot, more gives higher resolution. Default `n_points = 1000`.
- `u_name`: A `Symbol` giving the name of the field corresponding to the x component of the velocity. Default `u_name = :u`.
- `v_name`: A `Symbol` giving the name of the field corresponding to the x component of the velocity. Default `v_name = :v`.
"""
function check_windwater(
    windwater_itp::InterpolatedField,
    time::Real = 0.0;
    limits::NTuple{4, Real} = (-100, -50, 5, 35),
    n_points::Integer = 1000,
    u_name::Symbol = :u,
    v_name::Symbol = :v)

    xs = range(start=limits[1], stop=limits[2], length = n_points)
    ys = range(start=limits[3], stop=limits[4], length = n_points)
    us = [windwater_itp.fields[u_name](sph2xy(x, y, windwater_itp.ref)..., time) for x in xs, y in ys]
    vs = [windwater_itp.fields[v_name](sph2xy(x, y, windwater_itp.ref)..., time) for x in xs, y in ys]

    fig = default_fig()

    ax_u = geo_axis(fig[1, 1], title = "u", limits = limits)
    heatmap!(ax_u, xs, ys, us)

    ax_v = geo_axis(fig[1, 2], title = "v", limits = limits)
    heatmap!(ax_v, xs, ys, vs)    

    return fig
end