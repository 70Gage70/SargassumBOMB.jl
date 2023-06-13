using GeoMakie, CairoMakie

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
x_range = range(start = x0 - 1, length = 5, stop = x0 + 1)
y_range = range(start = y0 - 1, length = 5, stop = y0 + 1)
clump_parameters = BOMParameters(ref)
spring_parameters = SpringParameters(k -> 20, 0.4)

@named RRaft = RectangularRaft(x_range, y_range, clump_parameters, spring_parameters)
RRaft = structural_simplify(RRaft)

t_range = (0.0, 150.0)
# t_range = (0.0, 100.0)

prob = ODEProblem(
    RRaft, 
    [],
    t_range, 
    []
)

@info "Solving model."

sol = solve(prob)

@variables COM(t)[1:2]
com = [sol[COM[1]] ;; sol[COM[2]]]
traj = xy2sph(com, ref)
lon_traj, lat_traj = (traj[:, 1], traj[:, 2]) 
times = sol.t

@info "Plotting results."

fig = Figure(resolution = (1920, 1080))

ga(fig, row, col, title) = GeoAxis(
    fig[row, col],
    dest = "+proj=eqc", # https://proj.org/en/9.2/operations/projections/eqc.html
    lonlims = (-100, -50),
    latlims = (5, 35),
    coastlines = true,
    title = title
)

ln = lines!(ga(fig, 1, 1, "traj"), lon_traj, lat_traj; color = times, linewidth = 4)
Colorbar(fig[1,2], ln, label = "Days")

fig