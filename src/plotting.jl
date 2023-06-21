using CairoMakie, GeoDatasets

const lsm = GeoDatasets.landseamask(; resolution = 'i', grid = 1.25);

function geo_axis(fig; fig_pos = [1, 1], limits = (-100, -50, 5, 35), title = "Test")
    Axis(
        fig[fig_pos[1], fig_pos[2]],
        limits = limits, 
        title = L"%$(title)",
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


function plot_traj(lon_traj, lat_traj, times)
    fig = Figure(
        resolution = (1920, 1080), 
        fontsize = 50,
        figure_padding = (5, 50, 5, 5));

    ax = geo_axis(fig, fig_pos = [1, 1], limits = (-100, -50, 5, 35))
    coastlines!(ax)

    ln = lines!(ax, lon_traj, lat_traj; color = times, linewidth = 4);
    Colorbar(
        fig[1,2], 
        ln, 
        label = L"\textrm{Days}",
        ticklabelsize = 40, 
        # ticks = [0, 50, 100, 150], 
        tickformat = values -> [L"%$(value)" for value in values],
        height = Relative(2/4)
    );

    return fig
end