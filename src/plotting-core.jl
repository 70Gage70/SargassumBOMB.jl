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
` of the data in `traj` scaled according to the sargassum content at week `week`.

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


# """
#     plot(bop::BOMBOptimizationProblem; high_accuracy, type, showprogress)

# Plot the [`BOMBOptimizationProblem`](@ref) in `bop` as well as the `SargassumFromAFAI` distributions for comparison at 
# the initial and final times provided by `bop.tspan`.

# ### Optional Arguments 

# The following arguments are passed directly to `simulate`, the result of which is used to make the plot.

# - `log_scale`: Whether to display the graphs using a logarithmic scale.
# - `show_coast`: Highlight the coastlines in each graph via [`coast!`](@ref). Default `false`.
# - `show_clouds`: Highlight clouds/missing data in each graph via [`clouds!`](@ref). Default `false`.
# """
# function plot(
#     bop::BOMBOptimizationProblem,
#     ymw1::NTuple{3, Integer},
#     ymw2::NTuple{3, Integer},
#     dists::Dict{Tuple{Int64, Int64}, SargassumFromAFAI.SargassumDistribution},
#     rtr_waterwind::RaftTrajectory;
#     log_scale::Bool = true,
#     show_coast::Bool = false,
#     show_clouds::Bool = false)

#     set_theme!(GEO_THEME())

#     fig = Figure(
#         # size = (1920, 1080), 
#         # size = (2420, 2320),
#         size = (2920, 2320),
#         fontsize = 50,
#         figure_padding = (5, 100, 5, 5))

#     limits = (-100, -40, 5, 35)

#     ### AFAI
#     # initial distribution 
#     dist_initial = dists[ymw1[1:2]]
#     ax = Axis(fig[1, 1], limits = limits, title = "AFAI initial $(monthname(ymw1[2])), week $(ymw1[3])")
#     SargassumFromAFAI.sarg!(ax, dist_initial, ymw1[3], log_scale = log_scale)
#     show_coast ? coast!(ax, dist_initial) : nothing
#     show_clouds ? clouds!(ax, dist_initial, ymw1[3]) : nothing
#     land!(ax)

#     # final distribution
#     dist_final = dists[ymw2[1:2]]
#     ax = Axis(fig[1, 2], limits = limits, title = "AFAI final $(monthname(ymw2[2])), week $(ymw2[3])")
#     SargassumFromAFAI.sarg!(ax, dist_final, ymw2[3], log_scale = log_scale)
#     show_coast ? coast!(ax, dist_final) : nothing
#     show_clouds ? clouds!(ax, dist_final, ymw2[3]) : nothing
#     land!(ax)

#     ### BOMB
#     rtr = bop.opt_rtr

#     # initial distribution 
#     tstart = ymw2time(ymw1...)
#     rtr_initial = time_slice(rtr, (tstart, tstart))
#     ax = Axis(fig[2, 1], limits = limits, title = "eBOMB initial [optim] $(monthname(ymw1[2])), week $(ymw1[3])")
#     trajectory_hist!(ax, rtr_initial, dist_initial, ymw1[3], log_scale = log_scale)
#     show_coast ? coast!(ax, dist_initial) : nothing
#     show_clouds ? clouds!(ax, dist_initial, ymw1[3]) : nothing
#     land!(ax)

#     # final distribution 
#     tspan_end = ymwspan2weekspan(ymw1, ymw2) |> x -> (ymw2time(x[end - 1]...), ymw2time(x[end]...))
#     rtr_final = time_slice(rtr, tspan_end)
#     ax = Axis(fig[2, 2], limits = limits, title = "eBOMB final [optim] $(monthname(ymw2[2])), week $(ymw2[3])")
#     trajectory_hist!(ax, rtr_final, dist_final, ymw2[3], log_scale = log_scale)
#     show_coast ? coast!(ax, dist_final) : nothing
#     show_clouds ? clouds!(ax, dist_final, ymw2[3]) : nothing
#     land!(ax)

#     ### WATERWIND COMPARISON

#     # initial distribution 
#     ax = Axis(fig[3, 1], limits = limits, title = "LEEWAY initial [optim] $(monthname(ymw1[2])), week $(ymw1[3])")
#     rtr_initial = time_slice(rtr_waterwind, (tstart, tstart))
#     trajectory_hist!(ax, rtr_initial, dist_initial, ymw1[3], log_scale = log_scale)
#     show_coast ? coast!(ax, dist_initial) : nothing
#     show_clouds ? clouds!(ax, dist_initial, ymw1[3]) : nothing
#     land!(ax)

#     # final distribution 
#     ax = Axis(fig[3, 2], limits = limits, title = "LEEWAY final [optim] $(monthname(ymw2[2])), week $(ymw2[3])")
#     rtr_final = time_slice(rtr_waterwind, tspan_end)
#     trajectory_hist!(ax, rtr_final, dist_final, ymw2[3], log_scale = log_scale)
#     show_coast ? coast!(ax, dist_final) : nothing
#     show_clouds ? clouds!(ax, dist_final, ymw2[3]) : nothing
#     land!(ax)

#     ### LABELS
#     ltx(x) = latexify(x, fmt = FancyNumberFormatter(4))

#     loss_ltx = ltx(bop.opt)
#     p_vals = [bop.params[param].optimizable ? bop.params[param].opt : bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES]
#     δ_opt, a_opt, σ_opt, A_spring_opt, λ_opt, μ_max_opt, m_opt, k_N_opt, T_min_opt, T_max_opt = ltx.(p_vals)

#     loss_ltx_comp = ltx(bop.loss_func.f(rtr_waterwind))
    
#     if bop.rhs == Leeway!
#         fig[-4,:] = Label(fig, L"\text{LEEWAY}")
#     elseif bop.rhs == Raft!
#         fig[-4,:] = Label(fig, L"\text{eBOMB}")
#     end

#     fig[-3,:] = Label(fig, L"[%$(bop.loss_func.name)] Loss(eBOMB) =  %$(loss_ltx), Loss(LEEWAY) =  %$(loss_ltx_comp)")
    
#     fig[-2,:] = Label(fig, L"Optimals (clumps): $\delta =$ %$(δ_opt), $\tau =$ %$(a_opt), $\sigma =$ %$(σ_opt)")
#     fig[-1,:] = Label(fig, L"Optimals (springs): $A_\text{spring} =$ %$(A_spring_opt), $λ =$ %$(λ_opt), $\mu_\text{max} =$ %$(μ_max_opt), $m =$ %$(m_opt), $k_N =$ %$(k_N_opt)")
#     fig[0,:] = Label(fig, L"Optimals (biology): $\mu_\text{max} =$ %$(μ_max_opt), $m =$ %$(m_opt), $k_N =$ %$(k_N_opt), $T_\text{min} =$ %$(T_min_opt), $T_\text{max} =$ %$(T_max_opt)")


#     return fig
# end