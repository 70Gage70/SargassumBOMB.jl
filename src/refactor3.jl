
using DifferentialEquations, JLD2

include("vector-field-files.jl")
include("parameters.jl")
include("plotting.jl")

##########################################################

ref = EquirectangularReference(lon0 = -75.0, lat0 = 10.0)


v_x, v_y, Dv_xDt, Dv_yDt, u_x, u_y, Du_xDt, Du_yDt, ω = load("itp.jld2", 
    "v_x", "v_y", "Dv_xDt", "Dv_yDt", "u_x", "u_y", "Du_xDt", "Du_yDt","ω")

# Computing the spring force
Fs(xy1, xy2, parameters) = spring_force(xy1, xy2, parameters)


function Clump!(du, u, p::ClumpParameters, t)
    x, y = u
    α, τ, R, f = (p.α, p.τ, p.R, p.f)
    du[1] = u_x(x, y, t) + τ * (
        R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t) + (f + R*ω(x, y, t)/3)*u_y(x, y, t)
    )
    du[2] = u_y(x, y, t) + τ * (
        R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t) - (f + R*ω(x, y, t)/3)*u_x(x, y, t)
    )
end

function Clump(u, p::ClumpParameters, t)
    x, y = u
    α, τ, R, f = (p.α, p.τ, p.R, p.f)

    return [
    u_x(x, y, t) + τ * (
        R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t) + (f + R*ω(x, y, t)/3)*u_y(x, y, t)
    ),
    u_y(x, y, t) + τ * (
        R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t) - (f + R*ω(x, y, t)/3)*u_x(x, y, t)
    )
    ]
end


function Raft!(du, u, p::RaftParameters, t)
    for i = 1:length(p.clumps)
        du[i, :] .= Clump(u[i, :], p.clumps[i], t)

        for j = 1:length(p.clumps)
            if j != i
                du[i, :] .+= Fs(u[i,:], u[j,:], p.springs[i, j])
            end
        end
    end
end

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

function plot_com(sol)
    com = [vec(sum(sol.u[t], dims = 1))/size(sol.u[1], 1) for t = 1:length(sol.t)]
    traj = xy2sph(com, ref)
    lon_traj, lat_traj = (traj[:, 1], traj[:, 2]) 
    times = sol.t

    return plot_traj(lon_traj, lat_traj, times)
end