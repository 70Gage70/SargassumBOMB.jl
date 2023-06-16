using DifferentialEquations, JLD2

include("vector-field-files.jl")
include("parameters.jl")
include("plotting.jl")

##########################################################

ref = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)

const water_itp = load("water_itp.jld", "water_itp")
const wind_itp = load("wind_itp.jld", "wind_itp")

# All the functions depending on wind and water vector fields.
# Note that `water_itp` and `wind_itp` must be loaded before using these.
v_x(x, y, t) = water_itp.u(x, y, t)
v_y(x, y, t) =  water_itp.v(x, y, t)
Dv_xDt(x, y, t) = MaterialDerivativeX(water_itp, x, y, t)
Dv_yDt(x, y, t) = MaterialDerivativeY(water_itp, x, y, t)
u_x(x, y, t, α) = (1 - α) * water_itp.u(x, y, t) + α * wind_itp.u(x, y, t)
u_y(x, y, t, α) = (1 - α) * water_itp.v(x, y, t) + α * wind_itp.v(x, y, t)
Du_xDt(x, y, t, α) = (1 - α) * MaterialDerivativeX(water_itp, x, y, t) + α * MaterialDerivativeX(wind_itp, x, y, t) 
Du_yDt(x, y, t, α) = (1 - α) * MaterialDerivativeY(water_itp, x, y, t) + α * MaterialDerivativeY(wind_itp, x, y, t) 
ω(x, y, t) = Vorticity(water_itp, x, y, t)

# v_x(x, y, t) = 1
# v_y(x, y, t) =  1
# Dv_xDt(x, y, t) = 1
# Dv_yDt(x, y, t) = 1
# u_x(x, y, t, α) = 1
# u_y(x, y, t, α) = 1
# Du_xDt(x, y, t, α) = 1 
# Du_yDt(x, y, t, α) = 1
# ω(x, y, t) = 1

# Computing the spring force
Fs(xy1, xy2, parameters) = spring_force(xy1, xy2, parameters)
# Fs(xy1, xy2, parameters) = [1.0, 1.0]


function Clump!(du, u, p::ClumpParameters, t)
    x, y = u
    α, τ, R, f = (p.α, p.τ, p.R, p.f)
    du[1] = u_x(x, y, t, α) + τ * (
        R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α)
    )
    du[2] = u_y(x, y, t, α) + τ * (
        R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α)
    )
end

function Clump(u, p::ClumpParameters, t)
    x, y = u
    α, τ, R, f = (p.α, p.τ, p.R, p.f)

    return [
    u_x(x, y, t, α) + τ * (
        R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α)
    ),
    u_y(x, y, t, α) + τ * (
        R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α)
    )
    ]
end


function Raft!(du, u, p::RaftParameters, t)
    for i = 1:length(p.clumps)
        du[i, :] .= Clump(u[i, :], p.clumps[i], t) + sum([Fs(u[i,:], u[j,:], p.springs[i, j]) for j = 1:length(p.clumps) if j != i])
    end
end

# function Raft!(du, u, p::RaftParameters, t)
#     for i = 1:length(p.clumps)
#         du[i, :] .= Clump(u[i, :], p.clumps[i], t)

#         for j = 1:length(p.clumps)
#             du[i, :] .+= Fs(u[i,:], u[j,:], p.springs[i, j])
#         end
#     end
# end

# function Raft!(du, u, p::RaftParameters, t)
#     for i = 1:length(p.clumps)
#         x = u[i, 1]
#         y = u[i, 2]

#         α, τ, R, f = (p.clumps[i].α, p.clumps[i].τ, p.clumps[i].R, p.clumps[i].f)

#         du[i, 1] = u_x(x, y, t, α) + τ * (
#             R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α)
#         )
#         du[i, 2] = u_y(x, y, t, α) + τ * (
#             R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α)
#         )

#         xyi = u[i, :]
#         for j = 1:length(p.clumps)
#             if j != i
#                 xyj = u[j, :]

#                 du[i, 1] += (p.springs[i,j]).k(norm(xyi - xyj))*((p.springs[i,j]).L/norm(xyi - xyj) - 1)*(xyi[1] - xyj[1])
#                 du[i, 2] += (p.springs[i,j]).k(norm(xyi - xyj))*((p.springs[i,j]).L/norm(xyi - xyj) - 1)*(xyi[2] - xyj[2])
#             end
#         end
#     end
# end

# xy0 = sph2xy(-64, 14, ref) 
# tspan = (0.0, 150.0)
# p = ClumpParameters(ref)
# prob = ODEProblem(Clump!, xy0, tspan, p)

# sol = solve(prob)
# traj = xy2sph(sol.u, ref)
# lon_traj, lat_traj = (traj[:, 1], traj[:, 2]) 
# times = sol.t

# plot_traj(lon_traj, lat_traj)

function get_prob(n, tmax)
    @info "Building Problem."
    x0, y0 = sph2xy(-64, 14, ref) 
    x_range = range(start = x0 - 5, length = n, stop = x0 + 5)
    y_range = range(start = y0 - 5, length = n, stop = y0 + 5)
    clump_parameters = ClumpParameters(ref)
    spring_parameters = SpringParameters(k -> 20, step(x_range))
    raft_parameters = RectangularRaftParameters(x_range, y_range, clump_parameters, spring_parameters)
    tspan_raft = (0.0, tmax)
    prob_raft = ODEProblem(Raft!, raft_parameters.xy0, tspan_raft, raft_parameters)

    return prob_raft
end

# @info "Solving problem."
# @time sol = solve(prob_raft)
# return sol

function plot_com(sol)
    com = [vec(sum(sol.u[t], dims = 1))/size(sol.u[1], 1) for t = 1:length(sol.t)]
    traj = xy2sph(com, ref)
    lon_traj, lat_traj = (traj[:, 1], traj[:, 2]) 
    times = sol.t

    return plot_traj(lon_traj, lat_traj, times)
end