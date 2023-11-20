"""
    default_fig()
"""
function default_fig()
    return Figure(
    resolution = (1920, 1080), 
    fontsize = 50,
    figure_padding = (5, 100, 5, 5));
end

"""
    geo_axis(fig_pos; title, limits, xticks, yticks)
"""
function geo_axis(
    fig_pos::GridPosition;
    title::AbstractString = L"\mathrm{Title}",
    limits::Tuple = (-100, -50, 5, 35),
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
    vector_field_t!(axis, vf, t; fieldnames)
"""
function vector_field_t!(
    axis::Axis,
    vf::InterpolatedField,
    t::Real;
    fieldnames::Tuple{Symbol, Symbol} = (:u, :v))

    lims = axis.limits.val
    f(x, y) = Point2f(
        vf.fields[fieldnames[1]](sph2xy(x, y, vf.ref)..., t), 
        vf.fields[fieldnames[2]](sph2xy(x, y, vf.ref)..., t))
    
    streamplot!(axis, f, lims[1]..lims[2], lims[3]..lims[4])
end

"""
    scalar_field_t!(axis, sf, t; fieldnames, n_points)
"""
function scalar_field_t!(
    axis::Axis,
    sf::InterpolatedField,
    t::Real;
    fieldname::Symbol = :u,
    n_points::Integer = 100)

    lims = axis.limits.val

    xs = range(start = lims[1], stop = lims[2], length = n_points)
    ys = range(start = lims[3], stop = lims[4], length = n_points)
    zs = [sf.fields[fieldname](sph2xy(x, y, sf.ref)..., t) for x in xs, y in ys]
    
    heatmap!(axis, xs, ys, zs)
end

"""
    trajectory_hist!(axis, traj, lon_bins, lat_bins; opts...)
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