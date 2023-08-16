include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

"""
    is_land(land_itp, x, y)

Return `true` or `false` according to whether the point `(x, y)` is on the land according to `land_itp`.

Assumes that `land_itp` has been constructed using the default interpolation (nearest neighbor) and hence simply checks if `land_itp.u(x, y) == 1.0`.
"""
function is_land(land_itp::StaticField2DInterpolantEQR, x::Real, y::Real)
    return land_itp.u(x, y) == 1.0
end


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
    zs = [is_land(land_itp, sph2xy(x, y, land_itp.ref)...) for x in xs, y in ys]

    fig = default_fig()
    ax = geo_axis(fig[1, 1], title="Land", limits=limits)
    heatmap!(ax, xs, ys, zs)

    return fig
end


"""
    die_land(land_itp)

Create a `DiscreteCallback` which kills clumps when they reach the shore and terminates the integration if no clumps remain.

### Arguments

- `land_itp`: A `StaticField2DInterpolantEQR` which gives the interpolated land locations.
"""
function die_land(land_itp::StaticField2DInterpolantEQR)

    function condition(u, t, integrator)
        return any([is_land(land_itp, u[2*i-1], u[2*i])  for i = 1:Integer(length(u)/2)])
    end

    function affect!(integrator)
        u = integrator.u
        t = integrator.t
        inds = findall([is_land(land_itp, u[2*i-1], u[2*i])  for i = 1:Integer(length(u)/2)])
        inds = [inds[i] - (i - 1) for i = 1:length(inds)]
        # if we have to delete multiple clumps in one step, deleting one clump will change the indices of the others.
        # since findall is sorted, after you delete the clump indexed by inds[1], then the clumps with indices >inds[1]
        # have their index decreased by 1, and so on

        if length(inds) == Integer(length(u)/2) # all clumps will be removed, so terminate
            terminate!(integrator)
        else
            for i in inds
                deleteat!(integrator, 2*i - 1) # e.g. index i = 2, delete the 3rd component (x coord of 2nd clump)
                deleteat!(integrator, 2*i - 1) # now the y coordinate is where the x coordinate was
                kill!(integrator.p, i, t) # remove the clump and relabel its connections
            end
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