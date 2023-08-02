using DifferentialEquations, ModelingToolkit
using JLD2
using ProgressMeter

include(joinpath(@__DIR__, "..", "parameters.jl"))
include(joinpath(@__DIR__, "..", "..", "interpolants", "interpolant-derivatives.jl"))

# itp_path = joinpath(@__DIR__, "..", "..", "interpolants", "ias")
# isdefined(@__MODULE__, :water_itp) || (const water_itp = load(joinpath(itp_path, "water_itp.jld2"), "water_itp"))
# isdefined(@__MODULE__, :wind_itp) || (const wind_itp = load(joinpath(itp_path, "wind_itp.jld2"), "wind_itp"))
# isdefined(@__MODULE__, :ref_itp) || (const ref_itp = water_itp.ref)

itp_path = joinpath(@__DIR__, "..", "..", "interpolants", "glorys")
isdefined(@__MODULE__, :water_itp) || (const water_itp = load(joinpath(itp_path, "water_itp.jld2"), "water_itp"))
isdefined(@__MODULE__, :wind_itp) || (const wind_itp = load(joinpath(itp_path, "wind_itp.jld2"), "wind_itp"))
isdefined(@__MODULE__, :temp_itp) || (const temp_itp = load(joinpath(itp_path, "temp_itp.jld2"), "temp_itp"))
isdefined(@__MODULE__, :no3_itp) || (const no3_itp = load(joinpath(itp_path, "no3_itp.jld2"), "no3_itp"))
isdefined(@__MODULE__, :ref_itp) || (const ref_itp = water_itp.ref)

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

# v(x, y, t) = [water_itp.u(x, y, t), water_itp.v(x, y, t)]
# Dv_Dt(x, y, t) = [MaterialDerivativeX(water_itp, x, y, t), MaterialDerivativeY(water_itp, x, y, t)]
# u(x, y, t) = [
#     (1 - α) * water_itp.u(x, y, t) + α * wind_itp.u(x, y, t), 
#     (1 - α) * water_itp.v(x, y, t) + α * wind_itp.v(x, y, t)]
# Du_Dt(x, y, t) = [
#     (1 - α) * MaterialDerivativeX(water_itp, x, y, t) + α * MaterialDerivativeX(wind_itp, x, y, t), 
#     (1 - α) * MaterialDerivativeY(water_itp, x, y, t) + α * MaterialDerivativeY(wind_itp, x, y, t)]
# ω(x, y, t) = Vorticity(water_itp, x, y, t)

# Computing the spring force
Fs(xy1, xy2, parameters) = spring_force(xy1, xy2, parameters)

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


# think of u[i, j, xy], i.e. u[3, 4, 1:2] as the [x y] coordinates of the clump in row 3, column 4
function Raft!(du, u, p::RaftParameters, t)
    α, τ, R, f = p.clumps.α, p.clumps.τ, p.clumps.R, p.clumps.f

    for i = 1:size(u, 1), j = 1:size(u, 2)
        x, y = u[i, j, :]

        du[i, j, 1] = u_x(x, y, t, α) + τ * (
            R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α)
        )
        du[i, j, 2] = u_y(x, y, t, α) + τ * (
            R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α)
        )

        du[i, j, :] += τ*sum(spring_force(u[i, j, :], u[m, n, :], p.springs) for (m, n) in p.connections[(i, j)])
    end
end

function RaftCOM(sol::AbstractArray{<:Real, 4})
    [
        [sum(u[:, :, 1])/length(u[:, :, 1]), sum(u[:, :, 2])/length(u[:, :, 2])] for u in sol.u
    ]
end

function Raftij(sol::AbstractArray{<:Real, 4}, i::Integer, j::Integer)
    [sol.u[k][i, j, 1:2] for k = 1:length(sol.u)]
end

### clump
# xy0 = sph2xy(-64, 14, ref_itp) 
# tspan = (0.0, 150.0)
# cp = ClumpParameters(ref_itp)
# clump_prob = ODEProblem(Clump!, xy0, tspan, cp)

# @time sol = solve(clump_prob)
# traj = xy2sph(sol.u, ref_itp)
# lon_traj, lat_traj = (traj[:, 1], traj[:, 2]) 
# times = sol.t

### raft
# x_range = range(-63, -62, step = 0.1)
# y_range = range(13, 15, step = 0.1)
# tspan = (0.0, 15.0)
# cp = ClumpParameters(ref_itp)
# spring_k = x -> 20
# rp = RaftParameters(x_range, y_range, cp, spring_k)
# raft_prob = ODEProblem(Raft!, rp.xy0, tspan, rp)

# @time sol = solve(raft_prob, saveat = 5.0)

# n_traj = 100

# # prob_fun tells the ensemble how to "start over", in this case the i'th run is the problem remade with the i'th initial condition and time range.
# function prob_func(prob, i, repeat)
#     remake(prob, u0 = rp.xy0 .+ rand(), tspan = tspan)
# end

# # output_func generally handles what is done with the solution at each iteration; here we only use it to keep track of progress
# progress = Progress(n_traj)

# function output_func(sol, i)
#     next!(progress)
#     sol, false
# end

# Eprob = EnsembleProblem(
#     raft_prob, 
#     prob_func = prob_func,
#     output_func = output_func
# )

# sim = solve(Eprob, EnsembleThreads(), trajectories = n_traj, saveat = 5.0);

# finish!(progress)


# plot_traj(lon_traj, lat_traj)

# function get_prob(n, tmax, ref = ref_itp)
#     @info "Building Problem."
#     x0, y0 = sph2xy(-64, 14, ref) 
#     x_range = range(start = x0 - 5, length = n, stop = x0 + 5)
#     y_range = range(start = y0 - 5, length = n, stop = y0 + 5)
#     clump_parameters = ClumpParameters(ref)
#     spring_parameters = SpringParameters(k -> 20, step(x_range))
#     raft_parameters = RectangularRaftParameters(x_range, y_range, clump_parameters, spring_parameters)
#     tspan_raft = (0.0, tmax)
#     prob_raft = ODEProblem(Raft!, raft_parameters.xy0, tspan_raft, raft_parameters)

#     return prob_raft
# end

# @info "Solving problem."
# @time sol = solve(prob_raft)
# return sol

# function plot_com(sol)
#     com = [vec(sum(sol.u[t], dims = 1))/size(sol.u[1], 1) for t = 1:length(sol.t)]
#     traj = xy2sph(com, ref)
#     lon_traj, lat_traj = (traj[:, 1], traj[:, 2]) 
#     times = sol.t

#     return plot_traj(lon_traj, lat_traj, times)
# end