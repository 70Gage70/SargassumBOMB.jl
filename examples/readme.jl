# installation
import Pkg
Pkg.Registry.add("General")
Pkg.Registry.add(Pkg.RegistrySpec(url = "https://github.com/70Gage70/SargassumRegistry"))
Pkg.add("Sargassum")

# basic usage
using Sargassum

ymw_initial = (2018, 3, 1)
ymw_final = (2018, 4, 1)
rtr = RaftParameters(ymw_initial, ymw_final) |> simulate

# plotting
using CairoMakie

dist_final = DIST_1718[(2018, 4)]

set_theme!(GEO_THEME())
fig = Figure()

ax = Axis(fig[1, 1], limits = (-100, -40, 5, 30), aspect = DataAspect(), title = "SargassumFromAFAI (satellite)")
sarg!(ax, dist_final, 1, log_scale = true)
land!(ax)

tspan_end = (ymw2time(ymw_final...) - 7, ymw2time(ymw_final...))
rtr_final = time_slice(rtr, tspan_end)
ax = Axis(fig[2, 1], limits = (-100, -40, 5, 30), aspect = DataAspect(), title = "SargassumBOMB (simulation)")
trajectory_hist!(ax, rtr_final, dist_final, 1, log_scale = true)
land!(ax)

fig