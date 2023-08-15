include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/statistic-methods.jl"))

"""
    is_shore(water_itp, x, y, t; tol)

Return `true` or `false` according to whether the point `(x, y, t)` is on the shore according to `water_itp`.

A point is considered to be on the shore if both `water_itp.u` and `water_itp.v` are less than `tol` at that point.
"""
function is_shore(
    water_itp::VectorField2DInterpolantEQR,
    x::Real,
    y::Real,
    t::Real;
    tol::Real = 0.1)

    return (abs(water_itp.u(x, y, t)) < tol) && (abs(water_itp.v(x, y, t)) < tol)
end

"""
    avoid_shore(water_itp; tol)

Construct a `DiscreteCallback` which terminates integration when both the `u` and `v` components of `water_itp` are less than `tol`.
"""
function avoid_shore(
    water_itp::VectorField2DInterpolantEQR;
    tol::Real=0.1,
    save_positions::NTuple{2,Bool}=(true, false))

    condition(u, t, integrator) = is_shore(water_itp, u..., t, tol = tol)
    affect!(integrator) = terminate!(integrator)
    return DiscreteCallback(condition, affect!, save_positions=save_positions)
end


"""
    check_shore(water_itp; t, tol, limits, n_points)

Construct a plot of the inferred shore locations from `water_itp`.

### Arguments

- `water_itp`: A `VectorField2DInterpolantEQR` which gives the interpolated currents.

### Optional Arguments

- `t`: The time at which the plot is made. Default `t = 0.0`.
- `tol`: If `water_itp.u` and `water_itp.v` are both less than `tol` in absolute value, then that point is considered to be shore. Default `tol = 0.1`.
- `limits`: The `(lon_min, lon_max, lat_min, lat_max)` boundaries of the plot. 
- `n_points`: The number of points to use in the plot, more gives higher resolution. Default `n_points = 1000`.
"""
function check_shore(
    water_itp::VectorField2DInterpolantEQR;
    t::Real=0.0,
    tol::Real=0.1,
    limits=(-100, -10, 5, 35),
    n_points=1000)

    xs = range(start=limits[1], stop=limits[2], length=n_points)
    ys = range(start=limits[3], stop=limits[4], length=n_points)
    zs = [is_shore(water_itp, sph2xy(x, y, water_itp.ref)..., t, tol = tol) for x in xs, y in ys]

    fig = default_fig()
    ax = geo_axis(fig[1, 1], title="Shore", limits=limits)
    heatmap!(ax, xs, ys, zs)

    return fig
end


"""
    die_shore(water_itp; tol)

Create a `DiscreteCallback` which kills clumps when they reach the shore.

### Arguments

- `water_itp`: A `VectorField2DInterpolantEQR` which gives the interpolated currents.

### Optional Arguments

- `tol`: If `water_itp.u` and `water_itp.v` are both less than `tol` in absolute value, then that point is considered to be shore. Default `tol = 0.1`.
"""
function die_shore(
    water_itp::VectorField2DInterpolantEQR;
    tol::Real=0.1)

    function condition(u, t, integrator)
        return any([is_shore(water_itp, u[2*i-1], u[2*i], t, tol = tol)  for i = 1:Integer(length(u)/2)])
    end

    function affect!(integrator)
        u = integrator.u
        t = integrator.t
        inds = findall([is_shore(water_itp, u[2*i-1], u[2*i], t, tol = tol)  for i = 1:Integer(length(u)/2)])
        inds = [inds[i] - (i - 1) for i = 1:length(inds)]
        # if we have to delete multiple clumps in one step, deleting one clump will change the indices of the others.
        # since findall is sorted, after you delete the clump indexed by inds[1], then the clumps with indices >inds[1]
        # have their index decreased by 1, and so on

        for i in inds
            deleteat!(integrator, 2*i - 1) # e.g. index i = 2, delete the 3rd component (x coord of 2nd clump)
            deleteat!(integrator, 2*i - 1) # now the y coordinate is where the x coordinate was
            kill!(integrator.p, i, t) # remove the clump and its connections; amounts to removing index i and subtracting 1 from all larger indices
        end

        println("indices $inds hit shore at time $t")
    end

    return DiscreteCallback(condition, affect!)
end

function grow_test(t_grow::Vector{<:Real})
    function condition(u, t, integrator)
        return any([abs(t - tg) < 1.0 for tg in t_grow])
    end

    function affect!(integrator)
        u = integrator.u
        t = integrator.t

        resize!(integrator, length(u) + 2)
        u[end - 1] = u[end - 1 - 2] + rand()
        u[end] = u[end - 2] + rand()

        grow!(integrator.p, t)

        println("growth")
    end

    return DiscreteCallback(condition, affect!)
end