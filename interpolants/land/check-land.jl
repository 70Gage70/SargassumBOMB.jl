include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/statistic-methods.jl"))

################################################################################

"""
    check_land(land_itp; t, tol, limits, n_points)

Construct a plot of the land locations from `land_itp`.

### Arguments

- `land_itp`: A `StaticField2DInterpolantEQR` which gives the interpolated location of the land.

### Optional Arguments

- `limits`: The `(lon_min, lon_max, lat_min, lat_max)` boundaries of the plot. Default `limits = (-180, 180, -90, 90)`.
- `n_points`: The number of points to use in each dimension of the plot, more gives higher resolution. Default `n_points = 1000`.
"""
function check_land(
    land_itp::StaticField2DInterpolantEQR;
    limits=(-180, 180, -90, 90),
    n_points=1000)

    xs = range(start=limits[1], stop=limits[2], length=n_points)
    ys = range(start=limits[3], stop=limits[4], length=n_points)
    zs = [land_itp.u(sph2xy(x, y, land_itp.ref)...) for x in xs, y in ys]

    fig = default_fig()
    ax = geo_axis(fig[1, 1], title="Land", limits=limits)
    heatmap!(ax, xs, ys, zs)

    return fig
end