using DifferentialEquations
using ModelingToolkit
using GeoMakie, CairoMakie

include("parameters.jl")
include("coordinates.jl")
include("vector-field-files.jl")
include("vector-field-methods.jl")

#############################################

# variables can have bounds, parameters can have distributions

@info "Loading interpolants."

ref = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

# construct_wind_itp(wind_file_default, ref)
# construct_water_itp(water_file_default, ref)

@load "water_itp.jld"
@load "wind_itp.jld" 

@info "Generating model."

@parameters α τ R f
@variables t, x(t), y(t)

v_x(x, y, t) = water_itp.u(x, y, t)
v_y(x, y, t) =  water_itp.v(x, y, t)
Dv_xDt(x, y, t) = MaterialDerivativeX(water_itp, x, y, t)
Dv_yDt(x, y, t) = MaterialDerivativeY(water_itp, x, y, t)
u_x(x, y, t, α) = (1 - α) * water_itp.u(x, y, t) + α * wind_itp.u(x, y, t)
u_y(x, y, t, α) = (1 - α) * water_itp.v(x, y, t) + α * wind_itp.v(x, y, t)
Du_xDt(x, y, t, α) = (1 - α) * MaterialDerivativeX(water_itp, x, y, t) + α * MaterialDerivativeX(wind_itp, x, y, t) 
Du_yDt(x, y, t, α) = (1 - α) * MaterialDerivativeY(water_itp, x, y, t) + α * MaterialDerivativeY(wind_itp, x, y, t) 
ω(x, y, t) = Vorticity(water_itp, x, y, t)

@register_symbolic v_x(x, y, t)
@register_symbolic v_y(x, y, t)
@register_symbolic Dv_xDt(x, y, t)
@register_symbolic Dv_yDt(x, y, t)
@register_symbolic u_x(x, y, t, α)
@register_symbolic u_y(x, y, t, α)
@register_symbolic Du_xDt(x, y, t, α)
@register_symbolic Du_yDt(x, y, t, α)
@register_symbolic ω(x, y, t)

ddt = Differential(t)

@named BOM1 = ODESystem([
    ddt(x) ~ u_x(x, y, t, α) + τ * (
        R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α)
    ),
    ddt(y) ~ u_y(x, y, t, α) + τ * (
        R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α)
    )
])

# ics = sph2xy(-79.5, 25.5, ref) # gulf stream, use t_range = (0.0, 10.0)
ics = sph2xy(-64, 14, ref) # loop current, use t_range = (0.0, 200.0)
initial_conditions = [x => ics[1], y => ics[2]]
t_range = (0.0, 200.0)
params = BOM_parameters()
params = [α => params[1], τ => params[2], R => params[3], f => params[4]]

prob = ODEProblem(
    BOM1, 
    initial_conditions, 
    t_range, 
    params
)

@info "Solving model."

sol = solve(prob)
traj = xy2sph(sol.u, ref)
lon_traj, lat_traj = (traj[:,1], traj[:, 2]) #
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


# wind_itp = VectorField2DGridSPH(wind_file_default, lon_alias = "Lon", lat_alias = "Lat", lon_lat_time_order = [2, 1, 3])
# wind_itp = VectorField2DInterpolantSPH(wind_itp)

# lon_wind, lat_wind, t_wind, u_wind, v_wind = (wind_itp.lon, wind_itp.lat, wind_itp.time, wind_itp.u, wind_itp.v)
# u_wind_itp = [u_wind(lon, lat, 1) for lon in lon_wind, lat in lat_wind]
# v_wind_itp = [v_wind(lon, lat, 1) for lon in lon_wind, lat in lat_wind]

# surface!(ga(fig, 2, 1, "wind_x"), lon_wind, lat_wind, u_wind_itp[:,:,1])
# surface!(ga(fig, 2, 2, "wind_y"), lon_wind, lat_wind, v_wind_itp[:,:,1])

# water_itp = VectorField2DGridSPH(water_file_default, lon_lat_time_order = [2, 1, 3])
# water_itp = VectorField2DInterpolantSPH(water_itp)

# lon_water, lat_water, t_water, u_water, v_water = (water_itp.lon, water_itp.lat, water_itp.time, water_itp.u, water_itp.v)
# u_water_itp = [u_water(lon, lat, 1) for lon in lon_water, lat in lat_water]
# v_water_itp = [v_water(lon, lat, 1) for lon in lon_water, lat in lat_water]

# surface!(ga(fig, 3, 1, "water_x"), lon_water, lat_water, u_water_itp[:,:,1])
# surface!(ga(fig, 3, 2, "water_y"), lon_water, lat_water, v_water_itp[:,:,1])

fig