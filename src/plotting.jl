using GLMakie # use CairoMakie for vector graphics (slower)
using GeoDatasets

include("vector-fields/vector-field-files.jl")

########################################################

const lsm = GeoDatasets.landseamask(; resolution = 'i', grid = 1.25);

function geo_axis(fig; fig_pos = [1, 1], limits = (-100, -50, 5, 35), title = "Test")
    Axis(
        fig[fig_pos[1], fig_pos[2]],
        limits = limits, 
        title = title,
        xticklabelsize = 40,
        yticklabelsize = 40,
        xtickformat = values -> [
            if value > 0 
                L"%$(abs(value))^\circ \, \mathrm{E}" 
            elseif value == 0 
                L"0^\circ"
            elseif value < 0
                L"%$(abs(value))^\circ \, \mathrm{W}" 
            end
        for value in values],
        ytickformat = values -> [
            if value > 0 
                L"%$(abs(value))^\circ \, \mathrm{N}" 
            elseif value == 0 
                L"0^\circ"
            elseif value < 0
                L"%$(abs(value))^\circ \, \mathrm{S}" 
            end
        for value in values]
    )
 end

function coastlines!(axis)
    lon, lat, data = lsm
    contour!(axis, lon, lat, data, levels = [0.5], color = :black)
end

function land!(axis)
    lon, lat, data = lsm
    heatmap!(axis, lon, lat, data, colorrange = (0.3, 0.4), lowclip = RGBAf(1.0,1.0,1.0,0.0), highclip = RGBAf(0.6,0.6,0.6,1.0))
end

function arrows_timeslice(ax, t::Real, vf::VectorField2DInterpolantEQR)
    u = [vf.u(x, y, t) for x in vf.x, y in vf.y]
    v = [vf.v(x, y, t) for x in vf.x, y in vf.y]
    x, y = xy2sph(vf.x, vf.y, vf.ref)

    strength = vec(sqrt.(u .^ 2 .+ v .^ 2)) # color is proportional to norm

    return arrows!(ax, x, y, u, v,
        lengthscale = 0.01, # color indicates strength, not length of actual arrow
        arrowcolor = strength, linecolor = strength)
end

function streamlines_timeslice(ax, t::Real, vf::VectorField2DInterpolantEQR)
    function streamfunc(lon, lat; t = t)
        x, y = sph2xy(lon, lat, vf.ref)
        return Point2f(vf.u(x, y, t), vf.v(x, y, t))
    end

    lon, lat = xy2sph(water_itp.x, water_itp.y, water_itp.ref)

    return streamplot!(ax, streamfunc, lon, lat,
        density = 10.0, 
        # arrowsize = 6.0, linewidth = 6.0, 
        colormap = :buda)
end

function trajectory(ax, lons, lats, times; linewidth = 4)
    return lines!(ax, lons, lats; color = times, linewidth = linewidth)
end

function trajectory(ax, sol, ref::EquirectangularReference; linewidth = 4)
    traj = xy2sph(sol.u, ref)
    lons, lats = (traj[:, 1], traj[:, 2])
    return trajectory(ax, lons, lats, sol.t, linewidth = linewidth) 
end

function colorbar(fig, line; fig_pos = [1, 2])    
    return Colorbar(
        fig[fig_pos[1], fig_pos[2]], 
        line, 
        label = L"\textrm{Days}",
        ticklabelsize = 40, 
        # ticks = [0, 50, 100, 150], 
        tickformat = values -> [L"%$(value)" for value in values],
        height = Relative(2/4)
    )
end   

function geo_fig(;title = "Test")
    fig = Figure(
        resolution = (1920, 1080), 
        fontsize = 50,
        figure_padding = (5, 50, 5, 5));

    ax = geo_axis(fig, fig_pos = [1, 1], limits = (-100, -50, 5, 35), title = title)
    land!(ax)
    # coastlines!(ax)
    
    
    return (fig, ax)
end

# function plot_traj(lon_traj, lat_traj, times)
#     fig = Figure(
#         resolution = (1920, 1080), 
#         fontsize = 50,
#         figure_padding = (5, 50, 5, 5));

#     ax = geo_axis(fig, fig_pos = [1, 1], limits = (-100, -50, 5, 35))
#     coastlines!(ax)

#     ln = lines!(ax, lon_traj, lat_traj; color = times, linewidth = 4);
#     Colorbar(
#         fig[1,2], 
#         ln, 
#         label = L"\textrm{Days}",
#         ticklabelsize = 40, 
#         # ticks = [0, 50, 100, 150], 
#         tickformat = values -> [L"%$(value)" for value in values],
#         height = Relative(2/4)
#     );

#     return fig
# end