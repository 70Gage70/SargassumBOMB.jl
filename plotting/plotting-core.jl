"""
    default_fig()

Create a `Makie.Figure` with a size of `(1920, 1080)`, a fontsize of `50` and a padding 
of `(5, 100, 5, 5)`.
"""
function default_fig()
    return Figure(
    size = (1920, 1080), 
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

### Optional Arguments

- `colormap`: The colormap used in the colorbar. Default `SargassumColors.EUREKA`.
- `label_fontsize`: The font size of the label. Default 40.
- `tick_fontsize`: The font size of the label ticks. Default 40.
- `ticks`: A `Vector` of tick values. Default `[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]`.
- `barlength`: The length of the colorbar as a `Relative` size of the grid. Default `Relative(9/10).`
- `barwidth`: The width of the colorbar as a `Relative` size of the grid. Default `Relative(3/10).`
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
        # colormap = Reverse(:RdYlGn),
        colormap = EUREKA,
        colorscale = log10)
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
        # colormap = Reverse(:RdYlGn),
        colormap = EUREKA,
        colorscale = log10)
    )

    trajectory_hist!(axis, [traj], lon_bins, lat_bins; opts = opts)
end

"""
    trajectory_hist!(axis, traj, dist, week; opts...)

Create a `Makie.heatmap` on `axis` with the same bins as the `SargassumFromAFAI.SargassumDistribution` in `dist` 
` of the data in `traj` scaled according to the sargassum content at week `week`.

`traj` can be a single [`RaftTrajectory`](@ref) or a `Vector` of [`RaftTrajectory`](@ref). In the case of 
a `Vector`, all trajectories are mixed together to make a single plot.
"""
function trajectory_hist!(
    axis::Axis, 
    traj::Vector{<:RaftTrajectory},
    dist::SargassumDistribution,
    week::Integer;
    opts::NamedTuple = (
        # colormap = Reverse(:RdYlGn),
        colormap = EUREKA,
        colorscale = log10)
    )

    @assert week in [1, 2, 3, 4]

    lon_centers = dist.lon
    lat_centers = dist.lat

    binned = zeros(typeof(traj[1]).parameters[2], length(lon_centers), length(lat_centers))
    for tr in traj
        binned = binned + bins(tr, dist)
    end

    # rescale the trajectory data to be on the same scale as the distribution
    sarg = dist.sargassum[:,:,week]
    binned = binned*sum(sarg)/sum(binned)
    sarg_limits = (minimum(filter(x -> x > 0, sarg)), maximum(sarg))

    range_opts = (lowclip = :white, colorrange = sarg_limits)

    return heatmap!(axis, 
        lon_centers, 
        lat_centers, 
        binned; merge(opts, range_opts)...)
end

function trajectory_hist!(
    axis::Axis, 
    traj::RaftTrajectory,
    dist::SargassumDistribution,
    week::Integer;
    opts::NamedTuple = (
        # colormap = Reverse(:RdYlGn),
        colormap = EUREKA,
        colorscale = log10)
    )

    trajectory_hist!(axis, [traj], dist, week; opts = opts) 
end

"""
    plot(bop::BOMBOptimizationProblem; high_accuracy, type, showprogress)

Plot the [`BOMBOptimizationProblem`](@ref) in `bop` as well as the `SargassumFromAFAI` distributions for comparison at 
the initial and final times provided by `bop.tspan`.

### Optional Arguments 

The following arguments are passed directly to `simulate`, the result of which is used to make the plot.

- `opt`: .
"""
function plot(
    bop::BOMBOptimizationProblem,
    ymw1::NTuple{3, Integer},
    ymw2::NTuple{3, Integer},
    dists::Dict{Tuple{Int64, Int64}, SargassumFromAFAI.SargassumDistribution};
    opt::Bool = false)

    fig = Figure(
        # size = (1920, 1080), 
        # size = (2420, 2320),
        size = (2920, 2320),
        fontsize = 50,
        figure_padding = (5, 100, 5, 5))

    limits = (-100, -40, 5, 35)

    ### AFAI
    # initial distribution 
    dist_initial = dists[ymw1[1:2]]
    ax = geo_axis(fig[1, 1], limits = limits, title = "AFAI initial $(monthname(ymw1[2])), week $(ymw1[3])")
    SargassumFromAFAI.plot!(ax, dist_initial, ymw1[3], log_scale = true)
    land!(ax)

    # final distribution
    dist_final = dists[ymw2[1:2]]
    ax = geo_axis(fig[1, 2], limits = limits, title = "AFAI final $(monthname(ymw2[2])), week $(ymw2[3])")
    SargassumFromAFAI.plot!(ax, dist_final, ymw2[3], log_scale = true)
    land!(ax)

    ### BOMB
    rtr = bop.opt_rtr

    # initial distribution 
    tstart = ymw2time(ymw1...)
    rtr_initial = time_slice(rtr, (tstart, tstart))
    ax = geo_axis(fig[2, 1], limits = limits, title = "BOMB initial [optim] $(monthname(ymw1[2])), week $(ymw1[3])")
    trajectory_hist!(ax, rtr_initial, dist_initial, ymw1[3])
    land!(ax)

    # final distribution 
    tspan_end = ymwspan2weekspan(ymw1, ymw2)[end]
    rtr_final = time_slice(rtr, tspan_end)
    ax = geo_axis(fig[2, 2], limits = limits, title = "BOMB final [optim] $(monthname(ymw2[2])), week $(ymw2[3])")
    trajectory_hist!(ax, rtr_final, dist_final, ymw2[3])
    land!(ax)

    ### COMPARISON
    # δ_param = OptimizationParameter("δ",                1.25,   (1.05, 1.5),        false)
    # a_param = OptimizationParameter("a",                1.0e-4, (1.0e-5, 1.0e-3),   false)
    # σ_param = OptimizationParameter("σ",                1.0,    (0.9, 1.1),         false)
    δ_param = OptimizationParameter("δ",                3.0,   (2.99, 3.01),        false)
    a_param = OptimizationParameter("a",                4.25e-4, (4.24e-4, 4.26e-4),   false)
    σ_param = OptimizationParameter("σ",                0.0,    (-0.01, 0.01),         false)

    A_spring_param = OptimizationParameter("A_spring",  1.0,    (0.1, 3.0),         false)
    λ_param = OptimizationParameter("λ",                1.0,    (0.5, 1.5),         false)
    μ_max_param = OptimizationParameter("μ_max",        0.1,    (0.05, 0.5),        false)
    m_param = OptimizationParameter("m",                0.05,   (0.0, 0.1),         false)
    k_N_param = OptimizationParameter("k_N",            0.012,  (0.005, 0.05),      false)
    
    params_waterwind = Dict(param.name => param for param in [δ_param, a_param, σ_param, A_spring_param, λ_param, μ_max_param, m_param, k_N_param])
    
    bop_waterwind = BOMBOptimizationProblem(
        params = params_waterwind,
        rhs = WaterWind!,
        immortal = true,
        tspan = bop.tspan,
        n_levels = bop.n_levels,
        t_extra = bop.t_extra,
        loss_func = bop.loss_func
    )

    rtr = simulate(bop_waterwind, high_accuracy = high_accuracy, type = "default", showprogress = false)

    # initial distribution 
    ax = geo_axis(fig[3, 1], limits = limits, title = "WATER+WIND 3% initial [optim] $(monthname(start_date[2])), week 1")
    rtr_initial = time_slice(rtr, (tstart, tstart))
    trajectory_hist!(ax, rtr_initial, dist_initial, 1)
    land!(ax)

    # final distribution 
    ax = geo_axis(fig[3, 2], limits = limits, title = "WATER+WIND 3% final [optim] $(monthname(end_date[2])), week 1")
    rtr_final = time_slice(rtr, (tend - bop.t_extra, tend))
    trajectory_hist!(ax, rtr_final, dist_final, 1)
    land!(ax)

    # strings
    ltx(x) = latexify(x, fmt = FancyNumberFormatter(4))

    if type == "opt"
        loss_ltx = latexify(bop.opt, fmt = FancyNumberFormatter(4))
        p_vals = [bop.params[param].optimizable ? bop.params[param].opt : bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES]
    elseif type == "default"
        loss_ltx = latexify(loss(rtr, bop), fmt = FancyNumberFormatter(4))
        p_vals = [bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES]
    end

    loss_ltx_comp = latexify(loss(rtr_final, bop), fmt = FancyNumberFormatter(4))

    δ_opt, a_opt, σ_opt, A_spring_opt, λ_opt, μ_max_opt, m_opt, k_N_opt = ltx.(p_vals)

    if bop.rhs == WaterWind!
        fig[-4,:] = Label(fig, L"\text{WaterWind}")
    elseif bop.rhs == Raft!
        fig[-4,:] = Label(fig, L"\text{BOMB}")
    end

    fig[-3,:] = Label(fig, L"[%$(bop.loss_func.name)] Loss(BOMB) =  %$(loss_ltx), Loss(WATER) =  %$(loss_ltx_comp)")
    
    fig[-2,:] = Label(fig, L"Optimals (clumps): $\delta =$ %$(δ_opt), $a =$ %$(a_opt), $\sigma =$ %$(σ_opt)")
    fig[-1,:] = Label(fig, L"Optimals (springs): $A_\text{spring} =$ %$(A_spring_opt), $λ =$ %$(λ_opt), $\mu_\text{max} =$ %$(μ_max_opt), $m =$ %$(m_opt), $k_N =$ %$(k_N_opt)")
    fig[0,:] = Label(fig, L"Optimals (biology): $\mu_\text{max} =$ %$(μ_max_opt), $m =$ %$(m_opt), $k_N =$ %$(k_N_opt)")


    return fig
end