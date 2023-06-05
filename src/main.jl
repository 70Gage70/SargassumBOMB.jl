using DifferentialEquations
using ModelingToolkit
using GeoMakie, CairoMakie

include("parameters.jl")
include("coordinates.jl")
include("vector-fields.jl")

#############################################

# variables can have bounds, parameters can have distributions

@parameters α τ R f
@variables t x(t) y(t)

vfs = vector_fields()
# v_x(x, y, t) = vfs[1](x, y, t)
# v_y(x, y, t) = vfs[2](x, y, t)
# Dv_xDt(x, y, t) = vfs[3](x, y, t)
# Dv_yDt(x, y, t) = vfs[4](x, y, t)
# u_x(x, y, t, α) = vfs[5](x, y, t, α)
# u_y(x, y, t, α) = vfs[6](x, y, t, α)
# Du_xDt(x, y, t, α) = vfs[7](x, y, t, α)
# Du_yDt(x, y, t, α) = vfs[8](x, y, t, α)
# ω(x, y, t) = vfs[9](x, y, t)

lon0 = -75
lat0 = 10.0
v_x(x, y, t) = vfs[1](xy2sph(x, y, lon0, lat0)..., t)
v_y(x, y, t) = vfs[2](xy2sph(x, y, lon0, lat0)..., t)
Dv_xDt(x, y, t) = vfs[3](xy2sph(x, y, lon0, lat0)..., t)
Dv_yDt(x, y, t) = vfs[4](xy2sph(x, y, lon0, lat0)..., t)
u_x(x, y, t, α) = vfs[5](xy2sph(x, y, lon0, lat0)..., t, α)
u_y(x, y, t, α) = vfs[6](xy2sph(x, y, lon0, lat0)..., t, α)
Du_xDt(x, y, t, α) = vfs[7](xy2sph(x, y, lon0, lat0)..., t, α)
Du_yDt(x, y, t, α) = vfs[8](xy2sph(x, y, lon0, lat0)..., t, α)
ω(x, y, t) = vfs[9](xy2sph(x, y, lon0, lat0)..., t)

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

# initial_conditions = [x => -63, y => 14.0] # runs into Dominican
# t_range = (0.0, 5.4)
initial_conditions = [x => -83.4, y => 24.6]
t_range = (0.0, 1.0)
params = BOM_parameters()
params = [α => params[1], τ => params[2], R => params[3], f => params[4]]

prob = ODEProblem(
    BOM1, 
    initial_conditions, 
    t_range, 
    params
)

sol = solve(prob)
traj = stack(sol.u)
x_traj, y_traj = (traj[1,:], traj[2, :])
times = sol.t


fig = Figure(resolution = (1920, 1080))

ga(fig, row, col, title) = GeoAxis(
    fig[row, col],
    dest = "+proj=eqc", # https://proj.org/en/9.2/operations/projections/eqc.html
    lonlims = (-100, -50),
    latlims = (5, 35),
    coastlines = true,
    title = title
)


lines!(ga(fig, 1, 1, "traj"), x_traj, y_traj; color = times, linewidth = 4)



lon_wind, lat_wind, t_wind, u_wind, v_wind, lon_wtr, lat_wtr, t_wtr, u_wtr, v_wtr = load_vector_fields()
u_wtr_itp = [v_x(lon, lat, 1) for lon in lon_wtr, lat in lat_wtr]
v_wtr_itp = [v_y(lon, lat, 1) for lon in lon_wtr, lat in lat_wtr]

surface!(ga(fig, 2, 1, "water_x"), lon_wtr, lat_wtr, u_wtr[:,:,1])
surface!(ga(fig, 2, 2, "water_y"), lon_wtr, lat_wtr, v_wtr[:,:,1])

surface!(ga(fig, 3, 1, "wind_x_itp"), lon_wtr, lat_wtr, u_wtr_itp)
surface!(ga(fig, 3, 2, "wind_y_itp"), lon_wtr, lat_wtr, v_wtr_itp)

fig