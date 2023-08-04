include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../../CustomMakie.jl/src/statistic-methods.jl"))

"""
    avoid_shore(water_itp; tol)

Construct a `DiscreteCallback` which terminates integration when both the `u` and `v` components of `water_itp` are less than `tol`.
"""
function avoid_shore(
    water_itp::VectorField2DInterpolantEQR;
    tol::Real=0.1,
    save_positions::NTuple{2,Bool}=(true, false))

    condition(u, t, integrator) = (abs(water_itp.u(u..., t)) < tol) && (abs(water_itp.v(u..., t)) < tol)
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

    function is_shore(x, y, t)
        cond_u = abs(water_itp.u(sph2xy(x, y, water_itp.ref)..., t)) < tol
        cond_v = abs(water_itp.v(sph2xy(x, y, water_itp.ref)..., t)) < tol

        return cond_u && cond_v ? 1 : 0
    end

    zs = [is_shore(x, y, t) for x in xs, y in ys]

    fig = default_fig()
    ax = geo_axis(fig[1, 1], title="Shore", limits=limits)
    heatmap!(ax, xs, ys, zs)

    return fig
end

"""
    die_shore(n_clumps, water_itp; tol)

Create a `VectorContinuousCallback` which kills clumps when they reach the shore.

### Arguments

- `n_clumps`: An `Integer` which gives the number of clumps in the raft.
- `water_itp`: A `VectorField2DInterpolantEQR` which gives the interpolated currents.

### Optional Arguments

- `tol`: If `water_itp.u` and `water_itp.v` are both less than `tol` in absolute value, then that point is considered to be shore. Default `tol = 0.1`.
"""
function die_shore(
    n_clumps::Integer,
    water_itp::VectorField2DInterpolantEQR;
    tol::Real=0.1)

    function condition(out, u, t, integrator)
        for i = 1:Integer(length(u)/2)
            out[i] = max(max(abs(water_itp.u(u[2*i-1], u[2*i], t)), abs(water_itp.v(u[2*i-1], u[2*i], t))) - tol, 0.0)
            if out[i] == 0
                println("condition detonate at $i")
            end
        end
    end

    function affect!(integrator, idx)
        println("l1 = ", length(integrator.u))
        println("idx = ", idx)
        deleteat!(integrator, 2 * idx - 1) # e.g. idx = 2, delete the 3rd component (x coord of 2nd clump)
        deleteat!(integrator, 2 * idx - 1) # now the y coordinate is where the x coordinate was
        println("l2 = ", length(integrator.u))
        kill!(integrator.p, idx) # remove the clump and its connections
    end

    return VectorContinuousCallback(condition, affect!, n_clumps)
end

function die_shore2(
    water_itp::VectorField2DInterpolantEQR;
    tol::Real=0.1)

    function condition(u, t, integrator)
        return any(
            [max(abs(water_itp.u(u[2*i-1], u[2*i], t)), abs(water_itp.v(u[2*i-1], u[2*i], t))) < tol 
            for i = 1:Integer(length(u)/2)])
    end

    function affect!(integrator)
        println("l1 = ", length(integrator.u))
        deleteat!(integrator, 1) # e.g. idx = 2, delete the 3rd component (x coord of 2nd clump)
        deleteat!(integrator, 1) # now the y coordinate is where the x coordinate was
        println("l2 = ", length(integrator.u))
        kill!(integrator.p, 1) # remove the clump and its connections
    end

    return DiscreteCallback(condition, affect!)
end