using GeoMakie, CairoMakie, GeoDatasets

include("models.jl")
include("vector-field-files.jl")

###########################################

@info "Loading interpolants."

ref = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

# construct_wind_itp(wind_file_default, ref)
# construct_water_itp(water_file_default, ref)

@load "water_itp.jld"
@load "wind_itp.jld" 

@info "Generating model."

x0, y0 = sph2xy(-64, 14, ref) 
x_range = range(start = x0 - 5, length = 25, stop = x0 + 5)
y_range = range(start = y0 - 5, length = 25, stop = y0 + 5)
clump_parameters = BOMParameters(ref)
spring_parameters = SpringParameters(k -> 20, 0.4)

@named RRaft = RectangularRaft(x_range, y_range, clump_parameters, spring_parameters)
RRaft = structural_simplify(RRaft)

# t_range = (0.0, 150.0)
t_range = (0.0, 1.0)

@info "Generating problem."

prob = ODEProblem(
    RRaft, 
    [],
    t_range, 
    [],
    jac = true,
    sparse = true
)


# @info "Solving model."

# sol = solve(prob)

# @variables COM(t)[1:2]
# com = [sol[COM[1]] ;; sol[COM[2]]]
# traj = xy2sph(com, ref)
# lon_traj, lat_traj = (traj[:, 1], traj[:, 2]) 
# times = sol.t

# @info "Plotting results."

# function x_labeler(x::Real)
#     if x > 0.0
#         return L"%$(x)^\circ\,\textrm{E}"
#     elseif x == 0.0
#         return "0^\circ"
#     else if x < 0.0 
#         return L"%$(x)^\circ\,\textrm{W}"
#     end
# end

# function y_labeler(y::Real)
#     if y > 0.0
#         return "$(y)^\circ\,\textrm{N}"
#     elseif y == 0.0
#         return "0^\circ"
#     else if y < 0.0 
#         return "$(y)^\circ\,\textrm{S}"
#     end
# end

# fig = Figure(
#     resolution = (1920, 1080), 
#     fontsize = 50,
#     figure_padding = (5, 50, 5, 5));

# ax = Axis(
#     fig[1, 1],
#     limits = (-100, -50, 5, 35), 
#     title = L"\textrm{COM of a } 5 \times 5 \textrm{ raft}",
#     xticklabelsize = 40,
#     yticklabelsize = 40,
#     xtickformat = values -> [
#         if value > 0 
#             L"%$(abs(value))^\circ \, \mathrm{E}" 
#         elseif value == 0 
#             L"0^\circ"
#         elseif value < 0
#             L"%$(abs(value))^\circ \, \mathrm{W}" 
#         end
#     for value in values],
#     ytickformat = values -> [
#         if value > 0 
#             L"%$(abs(value))^\circ \, \mathrm{N}" 
#         elseif value == 0 
#             L"0^\circ"
#         elseif value < 0
#             L"%$(abs(value))^\circ \, \mathrm{S}" 
#         end
#     for value in values]
# );

# lon, lat, data = GeoDatasets.landseamask(; resolution = 'i', grid = 1.25);
# contour!(ax, lon, lat, data, levels = [0.5], color = :black);

# ln = lines!(ax, lon_traj, lat_traj; color = times, linewidth = 4);
# Colorbar(
#     fig[1,2], 
#     ln, 
#     label = L"\textrm{Days}",
#     ticklabelsize = 40, 
#     ticks = [0, 50, 100, 150], 
#     tickformat = values -> [L"%$(value)" for value in values],
#     height = Relative(2/4)
# );

# fig
