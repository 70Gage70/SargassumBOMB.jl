include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/statistic-methods.jl"))

################################################################################

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
        us = [itp.fields[u_name](sph2xy(x, y, itp.ref)..., time) for x in xs, y in ys]
        vs = [itp.fields[v_name](sph2xy(x, y, itp.ref)..., time) for x in xs, y in ys]

        ax_u = geo_axis(fig[i, 1], title = "$(names[i]) u", limits = limits)
        heatmap!(ax_u, xs, ys, us)

        ax_v = geo_axis(fig[i, 2], title = "$(names[i]) v", limits = limits)
        heatmap!(ax_v, xs, ys, vs)
    end    

    return fig
end