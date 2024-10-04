"""
    trajectory!(axis, xy, t; args...)

Add a `Makie.lines` plot to `axis` from the points `xy` and times `t`. Returns `Makie.lines!`.

`xy` can be a `Vector` of length `n_points` with entries that are vectors of length `2` or an `N x 2` matrix.

### Optional Arguments

- `args...`: All keyword arguments are passed directly to `Makie.lines!`.
"""
function trajectory!(
    axis::Axis, 
    xy::Vector{<:Vector{<:Real}},
    t::Vector{<:Real};
    args...)

    x, y = stack(xy, dims = 1)[:,1], stack(xy, dims = 1)[:,2]

    defaults = (color = t, linewidth = 2)

    return lines!(axis, x, y; merge(defaults, args)...)
end

function trajectory!(
    axis::Axis, 
    xy::Matrix{<:Real},
    t::Vector{<:Real};
    args...)

    x, y = xy[:,1], xy[:,2]

    defaults = (color = t, linewidth = 2)

    return lines!(axis, x, y; merge(defaults, args)...)
end


"""
    trajectory!(axis, traj; args...)

Add a `Makie.lines` plot to `axis` from the [`Trajectory`](@ref) in `traj`. Returns `Makie.lines!`.

### Optional Arguments

- `args...`: All keyword arguments are passed directly to `Makie.lines!`.
"""
function trajectory!(
    axis::Axis, 
    traj::Trajectory;
    args...)

    x, y = traj.xy[:,1], traj.xy[:,2]

    defaults = (color = traj.t, linewidth = 2)

    return lines!(axis, x, y; merge(defaults, args)...)
end


"""
    trajectory!(axis, rtraj; args...)

Add a `Makie.lines` plot to `axis` for each [`Trajectory`](@ref) in the [`RaftTrajectory`](@ref) in `rtraj`. Returns `nothing`.

### Optional Arguments

- `args...`: All keyword arguments are passed directly to `Makie.lines!`.
"""
function trajectory!(
    axis::Axis, 
    rtraj::RaftTrajectory;
    args...
    )

    for i in keys(rtraj.trajectories)
        x, y = rtraj.trajectories[i].xy[:,1], rtraj.trajectories[i].xy[:,2]

        defaults = (color = rtraj.trajectories[i].t, linewidth = 2)
        lines!(axis, x, y; merge(defaults, args)...)
    end
    
    return nothing
end

"""
    trajectory(rtr::RaftTrajectory; limits = (-100, -40, 5, 35))

Visualize a [`RaftTrajectory`](@ref) quickly.
"""
function trajectory(rtraj::RaftTrajectory; limits::NTuple{4, Real} = (-100, -40, 5, 35))
    set_theme!(GEO_THEME())
    fig = Figure()
    ax = Axis(fig[1, 1], limits = limits)

    for i in keys(rtraj.trajectories)
        x, y = rtraj.trajectories[i].xy[:,1], rtraj.trajectories[i].xy[:,2]

        lines!(ax, x, y, color = rtraj.trajectories[i].t, linewidth = 2)
    end

    limits = rtraj.com.t |> x -> x .- first(x) |> extrema
    t0 = time2datetime(rtraj.com.t[1])
    Colorbar(fig[1, 2], limits = limits, label = "Days since $(t0)")

    land!(ax)

    return fig
end

"""
    trajectory_hist!(axis, traj, lon_bins, lat_bins; log_scale, args...)

Create a `Makie.heatmap` on `axis` with bin centers at the coordinates defined by `lon_bins` 
and `lat_bins` of the data in `traj`.

`traj` can be a single [`RaftTrajectory`](@ref) or a `Vector` of [`RaftTrajectory`](@ref). In the case of 
a `Vector`, all trajectories are mixed together to make a single plot.

Returns `Makie.heatmap!`.

### Optional Arguments

- `log_scale`: Plot on a `log10` scale. Default `false`.
- `args...`: All keyword arguments are passed directly to `Makie.heatmap!`.
"""
function trajectory_hist!(
    axis::Axis, 
    traj::Vector{<:RaftTrajectory},
    lon_bins::StepRangeLen,
    lat_bins::StepRangeLen;
    log_scale::Bool = false,
    args...
    )

    δ_lon = step(lon_bins)
    δ_lat = step(lat_bins)

    lon_centers = [lon + δ_lon/2 for lon in collect(lon_bins)[1:end-1]]
    lat_centers = [lat + δ_lat/2 for lat in collect(lat_bins)[1:end-1]]

    binned = zeros(typeof(traj[1]).parameters[2], length(lon_bins) - 1, length(lat_bins) - 1)
    for tr in traj
        binned = binned + bins(tr, lon_bins, lat_bins)
    end

    defaults = (
        colormap = EUREKA,
        colorscale = log_scale ? log10 : x -> x,
        lowclip = :white,
        colorrange = (1.0, maximum(binned))
    )

    return heatmap!(axis, 
        lon_centers, 
        lat_centers, 
        binned; merge(defaults, args)...)
end

function trajectory_hist!(
    axis::Axis, 
    traj::RaftTrajectory,
    lon_bins::StepRangeLen,
    lat_bins::StepRangeLen;
    log_scale::Bool = false,
    args...
    )

    trajectory_hist!(axis, [traj], lon_bins, lat_bins; log_scale = log_scale, args...)
end

"""
    trajectory_hist!(axis, traj, dist, week; log_scale, args...)

Create a `Makie.heatmap` on `axis` with the same bins as the `SargassumFromAFAI.SargassumDistribution` in `dist` 
of the data in `traj` scaled according to the sargassum content at week `week`.

`traj` is a [`RaftTrajectory`](@ref).

Returns `Makie.heatmap!`.

### Optional Arguments

- `log_scale`: Plot on a `log10` scale. Default `false`.
- `args...`: All keyword arguments are passed directly to `Makie.heatmap!`.
"""
function trajectory_hist!(
    axis::Axis, 
    traj::RaftTrajectory,
    dist::SargassumDistribution,
    week::Integer;
    log_scale::Bool = false,
    args...
    )

    @argcheck week in [1, 2, 3, 4]

    lon_centers = dist.lon
    lat_centers = dist.lat
    binned = bins(traj, dist)

    # rescale the trajectory data to be on the same scale as the distribution
    sarg = dist.sargassum[:,:,week]
    binned = binned*sum(sarg)/sum(binned)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg))

    defaults = (
        colormap = EUREKA,
        colorscale = log_scale ? log10 : x -> x,
        lowclip = :white,
        colorrange = sarg_limits
    )

    return heatmap!(axis, 
        lon_centers, 
        lat_centers, 
        binned; merge(defaults, args)...)
end