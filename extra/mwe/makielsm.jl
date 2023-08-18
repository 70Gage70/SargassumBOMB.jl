# using GeoDatasets, CairoMakie

# lon, lat, data = GeoDatasets.landseamask(; resolution = 'i', grid = 1.25)

# fig = Figure()
# ax = Axis(fig[1,1], limits = (-100, -50, 5, 35))

# contour!(ax, lon, lat, data, levels = [0.5])

# fig

include("../src/plotting.jl")

fig, ax = geo_fig(title = L"\mathrm{BOM} \, 100 \, \mathrm{day} \,  \mathrm{trajectories}")
fig

# @time fig
