"""
    default_fig()

Create a `Makie.Figure` with a resolution of `(1920, 1080)`, a fontsize of `50` and a padding 
of `(5, 100, 5, 5)`.
"""
function default_fig()
    return Figure(
    resolution = (1920, 1080), 
    fontsize = 50,
    figure_padding = (5, 100, 5, 5));
end

"""
    geo_axis(fig_pos; title, limits, xticks, yticks)

Create a `Makie.Axis` suitable for plotting on an equirectangular projection.

### Arguments

- `fig_pos`: A `Makie.GridPosition` where the plot should go. For example if `fig` is a `Makie.Figure`, \
then `fig_pos[1, 1]` puts the axis in the first row and first column of `fig`.

### Optional Arguments 

- `title`: An `AbstractString`. Default `L"\\mathrm{Title}"`, where `Makie.L"..."` creates a `LaTeXString`.
- `limits`: `An NTuple{4, <:Real}` of the form `(xmin, xmax, ymin, ymax)`.
- `xticks`: A list of x tick mark locations.
- `yticks`: A list of y tick mark locations.
"""
function geo_axis(
    fig_pos::GridPosition;
    title::AbstractString = L"\mathrm{Title}",
    limits::NTuple{4, <:Real} = (-100, -50, 5, 35),
    xticks::Vector{<:Real} = [limits[1], limits[2]],
    yticks::Vector{<:Real} = [limits[3], limits[4]]
    )

    Axis(
        fig_pos,
        limits = limits, 
        title = title,
        xticks = xticks,
        yticks = yticks,
        xticklabelsize = 40,
        yticklabelsize = 40,
        xtickformat = values -> [
            if value > 0 
                L"%$(abs(value)) \, \degree \mathrm{E}" 
            elseif value == 0 
                L"0\degree"
            elseif value < 0
                L"%$(abs(value)) \, \degree \mathrm{W}" 
            end
        for value in values],
        ytickformat = values -> [
            if value > 0 
                L"%$(abs(value)) \, \degree \mathrm{N}" 
            elseif value == 0 
                L"0\degree"
            elseif value < 0
                L"%$(abs(value)) \, \degree \mathrm{S}" 
            end
        for value in values]
    )
end

"""
    land!(axis; landpath, color)

Add a land polygon to `axis::Makie.Axis`. This will be placed on top of any graphics that 
are already on the axis.

### Optional Arguments

- `landpath`: A `String` pointing to the location of the `.geojson` file containing the land features. By \
default this is provided by a Natural Earth 50 m file \
"https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne\\_50m\\_land.geojson" \
which is included automatically. 
- `color`: A `Makie.Colorant` which gives the color of the land. Default is grey via `RGBf(0.5,0.5,0.5)`.
"""
function land!(
    axis::Axis; 
    landpath::String = joinpath(@__DIR__, "geojson/ne_50m_land.geojson"), 
    color::Colorant = RGBf(0.5,0.5,0.5)) # aka colorant"gray50"

    landpoly = GeoJSON.read(read(landpath, String))
    poly!(axis, landpoly, color = color)
end

"""
    data_legend!(fig_pos, label; colormap, label_fontsize, tick_fontsize, ticks, barlength, barwidth)

Add a `Makie.Colorbar` with label `label` to the `GridPosion` in `fig_pos.`
"""
function data_legend!(
    fig_pos::GridPosition,
    label::AbstractString = L"\mathrm{Label}";
    colormap::Union{Symbol, Reverse{Symbol}} = :viridis,
    label_fontsize::Real = 40,
    tick_fontsize::Real = 40,
    ticks::Vector{<:Real} = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    barlength::Relative = Relative(8/10),
    barwidth::Relative = Relative(3/20))

    data_legend = GridLayout(fig_pos)

    Label(
        data_legend[1,1], 
        label,
        # L"\left(\frac{\pi_{\text{stat}}}{%$(latexify(colors_max, fmt = FancyNumberFormatter(3), env = :raw))}\right)^{1/4}", 
        fontsize = label_fontsize,
        valign = :bottom, 
        tellheight = false
    )
    
    Colorbar(
            data_legend[2, 1], 
            limits = extrema(ticks),
            colormap = colormap,
            ticklabelsize = tick_fontsize, 
            ticks = ticks, 
            tickformat = values -> [latexify(value, fmt = FancyNumberFormatter(2)) for value in values],
            width = barwidth,
            valign = :top,
            tellheight = false
    )

    rowsize!(data_legend, 2, barlength) # relative height of bar and its label
    # colsize!(fig.layout, 2, Relative(1/10)) # relative width of data legend    
end

"""
    trajectory!(axis, xy, t; opts...)

Add a `Makie.lines` plot to `axis` from the points `xy` and times `t`. 

`xy` can be a `Vector` of length `n_points` with entries that are vectors of length `2` or an `N x 2` matrix.
"""
function trajectory!(
    axis::Axis, 
    xy::Vector{<:Vector{<:Real}},
    t::Vector{<:Real};
    opts::NamedTuple = (
        color = t,
        linewidth = 2)
    )

    x, y = stack(xy, dims = 1)[:,1], stack(xy, dims = 1)[:,2]

    return lines!(axis, x, y; opts...)
end

function trajectory!(
    axis::Axis, 
    xy::Matrix{<:Real},
    t::Vector{<:Real};
    opts::NamedTuple = (
        color = t,
        linewidth = 2)
    )

    x, y = xy[:,1], xy[:,2]

    return lines!(axis, x, y; opts...)
end


"""
    trajectory!(axis, traj; opts...)

Add a `Makie.lines` plot to `axis` from the [`Trajectory`](@ref) in `traj`.
"""
function trajectory!(
    axis::Axis, 
    traj::Trajectory;
    opts::NamedTuple = (
        color = traj.t,
        linewidth = 2)
    )

    x, y = traj.xy[:,1], traj.xy[:,2]

    return lines!(axis, x, y; opts...)
end

function trajectory!(
    axis::Axis, 
    rtraj::RaftTrajectory;
    opts::NamedTuple = (
        linewidth = 2,)
    )

    for i in keys(rtraj.trajectories)
        x, y = rtraj.trajectories[i].xy[:,1], rtraj.trajectories[i].xy[:,2]
        lines!(axis, x, y, color = rtraj.trajectories[i].t; opts...)
    end
    
    return nothing
end

"""
    trajectory_hist!(axis, traj, lon_bins, lat_bins; opts...)

Create a `Makie.heatmap` on `axis` with bin centers at the coordinates defined by `lon_bins` 
and `lat_bins` of the data in `traj`.

`traj` can be a single [`RaftTrajectory`](@ref) or a `Vector` of [`RaftTrajectory`](@ref). In the case of 
a `Vector`, all trajectories are mixed together to make a single plot.
"""
function trajectory_hist!(
    axis::Axis, 
    traj::Vector{<:RaftTrajectory},
    lon_bins::StepRangeLen,
    lat_bins::StepRangeLen;
    opts::NamedTuple = (
        colormap = Reverse(:RdYlGn),
        colorscale = x -> x == 0.0 ? -1.0 : log10(x))
    )

    δ_lon = step(lon_bins)
    δ_lat = step(lat_bins)

    lon_centers = [lon + δ_lon/2 for lon in collect(lon_bins)[1:end-1]]
    lat_centers = [lat + δ_lat/2 for lat in collect(lat_bins)[1:end-1]]

    binned = zeros(typeof(traj[1]).parameters[2], length(lon_bins) - 1, length(lat_bins) - 1)
    for tr in traj
        binned = binned + bins(tr, lon_bins, lat_bins)
    end

    range_opts = (lowclip = :white, colorrange = (1.0, maximum(binned)))

    return heatmap!(axis, 
        lon_centers, 
        lat_centers, 
        binned; merge(opts, range_opts)...)
end

function trajectory_hist!(
    axis::Axis, 
    traj::RaftTrajectory,
    lon_bins::StepRangeLen,
    lat_bins::StepRangeLen;
    opts::NamedTuple = (
        colormap = Reverse(:RdYlGn),
        colorscale = x -> x == 0.0 ? -1.0 : log10(x))
    )

    trajectory_hist!(axis, [traj], lon_bins, lat_bins; opts = opts)
end

"""
    trajectory_hist!(axis, traj, dist; opts...)

Create a `Makie.heatmap` on `axis` with the same bins as the `SargassumFromAFAI.SargassumDistribution` in `dist` 
` of the data in `traj`.

`traj` can be a single [`RaftTrajectory`](@ref) or a `Vector` of [`RaftTrajectory`](@ref). In the case of 
a `Vector`, all trajectories are mixed together to make a single plot.
"""
function trajectory_hist!(
    axis::Axis, 
    traj::Vector{<:RaftTrajectory},
    dist::SargassumDistribution;
    opts::NamedTuple = (
        colormap = Reverse(:RdYlGn),
        colorscale = x -> x == 0.0 ? -1.0 : x)
    )

    lon_centers = dist.lon
    lat_centers = dist.lat

    binned = zeros(typeof(traj[1]).parameters[2], length(lon_centers), length(lat_centers))
    for tr in traj
        binned = binned + bins(tr, dist)
    end

    range_opts = (lowclip = :white, colorrange = (1.0, maximum(binned)))

    return heatmap!(axis, 
        lon_centers, 
        lat_centers, 
        binned; merge(opts, range_opts)...)
end

function trajectory_hist!(
    axis::Axis, 
    traj::RaftTrajectory,
    dist::SargassumDistribution;
    opts::NamedTuple = (
        colormap = Reverse(:RdYlGn),
        colorscale = x -> x == 0.0 ? -1.0 : x)
    )

    trajectory_hist!(axis, [traj], dist; opts = opts) 
end

"""
    plot(bop::BOMBOptimizationProblem)
"""
function plot_bop(bop::BOMBOptimizationProblem)
    if bop.opt === nothing
        @warn "BOMBOptimizationProblem is not optimized; showing defaults."
    end

    initial_time = first(bop.tspan)
    final_time = last(bop.tspan)

    fig = Figure(
        # resolution = (1920, 1080), 
        resolution = (2420, 2320),
        fontsize = 50,
        figure_padding = (5, 100, 5, 5))

    limits = (-100, -40, 5, 35)

    ### AFAI
    # initial distribution (AFAI)
    ax = geo_axis(fig[1, 1], limits = limits, title = "AFAI initial $(monthname(initial_time[2])), week 1")
    SFA_plot!(ax, initial_time, 1)
    land!(ax)

    # final distribution (AFAI)
    ax = geo_axis(fig[1, 2], limits = limits, title = "AFAI final $(monthname(final_time[2])), week 1")
    SFA_plot!(ax, final_time, 1)
    land!(ax)

    ### UNOPTIMIZED
    # initial distribution (SIMUL, unoptimized)
    ax = geo_axis(fig[2, 1], limits = limits, title = "SIMUL initial [default] $(monthname(initial_time[2])), week 1")
    rtr_dt, tstart, tend = integrate_bomb(bop, type = "default")
    dist = DISTS_2018[initial_time]
    rtr_dt_initial = time_slice(rtr_dt, (first(rtr_dt.t), first(rtr_dt.t)))
    trajectory_hist!(ax, rtr_dt_initial, dist)
    land!(ax)

    # final distribution (SIMUL, unoptimized)
    ax = geo_axis(fig[2, 2], limits = limits, title = "SIMUL final [default] $(monthname(final_time[2])), week 1")
    rtr_final = time_slice(rtr_dt, (tend - 8, tend))
    trajectory_hist!(ax, rtr_final, dist)
    land!(ax)

    ### OPTIMIZED
    # initial distribution (SIMUL, unoptimized)
    ax = geo_axis(fig[3, 1], limits = limits, title = "SIMUL initial [optim] $(monthname(initial_time[2])), week 1")
    rtr_dt, tstart, tend = integrate_bomb(bop, type = "opt")
    dist = DISTS_2018[initial_time]
    rtr_dt_initial = time_slice(rtr_dt, (first(rtr_dt.t), first(rtr_dt.t)))
    trajectory_hist!(ax, rtr_dt_initial, dist)
    land!(ax)

    # final distribution (SIMUL, unoptimized)
    ax = geo_axis(fig[3, 2], limits = limits, title = "SIMUL final [optim] $(monthname(final_time[2])), week 1")
    rtr_final = time_slice(rtr_dt, (tend - 8, tend))
    trajectory_hist!(ax, rtr_final, dist)
    land!(ax)

    # strings
    default_loss = loss_bomb(bop, "default")
    optimized_loss = loss_bomb(bop, "opt")

    ltx(x) = latexify(x, fmt = FancyNumberFormatter(4))

    dl_ltx, ol_ltx = ltx(default_loss), ltx(optimized_loss)
    ol_ltx = latexify(optimized_loss, fmt = FancyNumberFormatter(4))
    δ_def, a_def, β_def, A_spring_def, μ_max_def, m_def, k_N_def = ltx.([bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES])
    δ_opt, a_opt, β_opt, A_spring_opt, μ_max_opt, m_opt, k_N_opt = ltx.([bop.params[param].opt for param in OPTIMIZATION_PARAMETER_NAMES])

    if bop.rhs == WaterWind!
        fig[-3,:] = Label(fig, L"\text{WaterWind}")
    elseif bop.rhs == Raft!
        fig[-3,:] = Label(fig, L"\text{BOMB}")
    end

    fig[-2,:] = Label(fig, L"[%$(bop.loss_func.name)] Loss(default) = %$(dl_ltx), Loss(opt) =  %$(ol_ltx)")

    fig[-1,:] = Label(fig, L"Defaults: $\delta =$ %$(δ_def), $a =$ %$(a_def), $\beta =$ %$(β_def), $A_\text{spring} =$ %$(A_spring_def), $\mu_\text{max} =$ %$(μ_max_def), $m =$ %$(m_def), $k_N =$ %$(k_N_def)")
    
    fig[0,:] = Label(fig, L"Optimals: $\delta =$ %$(δ_opt), $a =$ %$(a_opt), $\beta =$ %$(β_opt), $A_\text{spring} =$ %$(A_spring_opt), $\mu_\text{max} =$ %$(μ_max_opt), $m =$ %$(m_opt), $k_N =$ %$(k_N_opt)")

    outfile = joinpath(@__DIR__, "..", "figures", "opt_test.png")
    rm(outfile, force = true)
    save(outfile, fig)

    return fig
end