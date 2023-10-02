using OrdinaryDiffEq
using JLD2

include(joinpath(@__DIR__, "parameters.jl"))
include(joinpath(@__DIR__, "..", "interpolants", "itp-core.jl"))
include(joinpath(@__DIR__, "..", "interpolants", "itp-derivatives.jl"))

########################################################################

# loading interpolants
itp_path = joinpath(@__DIR__, "..", "interpolants", "ocean-atmos")
isdefined(@__MODULE__, :water_itp) || (const water_itp = load(joinpath(itp_path, "water_itp.jld2"), "water_itp"))
isdefined(@__MODULE__, :wind_itp) || (const wind_itp = load(joinpath(itp_path, "wind_itp.jld2"), "wind_itp"))

# All the functions depending on wind and water vector fields.
# v_x(x, y, t) = water_itp.fields[:u](x, y, t)
# v_y(x, y, t) =  water_itp.fields[:v](x, y, t)
# Dv_xDt(x, y, t) = MaterialDerivativeX(water_itp, x, y, t)
# Dv_yDt(x, y, t) = MaterialDerivativeY(water_itp, x, y, t)
# u_x(x, y, t, α) = (1 - α) * water_itp.fields[:u](x, y, t) + α * wind_itp.fields[:u](x, y, t)
# u_y(x, y, t, α) = (1 - α) * water_itp.fields[:v](x, y, t) + α * wind_itp.fields[:v](x, y, t)
# Du_xDt(x, y, t, α) = (1 - α) * MaterialDerivativeX(water_itp, x, y, t) + α * MaterialDerivativeX(wind_itp, x, y, t) 
# Du_yDt(x, y, t, α) = (1 - α) * MaterialDerivativeY(water_itp, x, y, t) + α * MaterialDerivativeY(wind_itp, x, y, t) 
# ω(x, y, t) = Vorticity(water_itp, x, y, t)

v_x(x, y, t) = water_itp.fields[:u](x, y, t)
v_y(x, y, t) =  water_itp.fields[:v](x, y, t)
Dv_xDt(x, y, t) = water_itp.fields[:DDt_x](x, y, t)
Dv_yDt(x, y, t) = water_itp.fields[:DDt_y](x, y, t)
u_x(x, y, t, α) = (1 - α) * water_itp.fields[:u](x, y, t) + α * wind_itp.fields[:u](x, y, t)
u_y(x, y, t, α) = (1 - α) * water_itp.fields[:v](x, y, t) + α * wind_itp.fields[:v](x, y, t)
Du_xDt(x, y, t, α) = (1 - α) * water_itp.fields[:DDt_x](x, y, t) + α * wind_itp.fields[:DDt_x](x, y, t) 
Du_yDt(x, y, t, α) = (1 - α) * water_itp.fields[:DDt_x](x, y, t) + α * wind_itp.fields[:DDt_y](x, y, t) 
ω(x, y, t) = water_itp.fields[:vorticity](x, y, t)

########################################################################
########################################################################
########################################################################

"""
    Raft!(du, u, p::RaftParameters, t)

Compute the right-hand-side of the differential equation controlling the motion of a raft with parameters 
given by [`RaftParameters`](@ref).

The solution vector `u` is a vector of length `2n_clumps + 1` such that `u[1]` is an "amount" parameter which 
controls the growth and death of clumps by biophysical effects. Then, `u[2*i:2*i+1]` for `i = 1:n_clumps` gives 
the `[x, y]` coordinates of the clump in position `i`.
"""
function Raft!(du, u, p::RaftParameters, t)
    du[1] = p.gd_model.dSdt(u, t)

    α, τ, R, f = p.clumps.α, p.clumps.τ, p.clumps.R, p.clumps.f

    for i = 1:floor(Int64, length(u)/2)
        # note that x, y = u[2*i], u[2*i+1]

        du[2*i] = 
            u_x(u[2*i], u[2*i+1], t, α) 
            + τ * (
            R*Dv_xDt(u[2*i], u[2*i+1], t) 
            - R*(f + ω(u[2*i], u[2*i+1], t)/3)*v_y(u[2*i], u[2*i+1], t) 
            - Du_xDt(u[2*i], u[2*i+1], t, α) 
            + (f + R*ω(u[2*i], u[2*i+1], t)/3)*u_y(u[2*i], u[2*i+1], t, α)
        )
        du[2*i+1] = 
            u_y(u[2*i], u[2*i+1], t, α) 
            + τ * (
            R*Dv_yDt(u[2*i], u[2*i+1], t) 
            + R*(f + ω(u[2*i], u[2*i+1], t)/3)*v_x(u[2*i], u[2*i+1], t) 
            - Du_yDt(u[2*i], u[2*i+1], t, α) 
            - (f + R*ω(u[2*i], u[2*i+1], t)/3)*u_x(u[2*i], u[2*i+1], t, α)
        )

        # @views combined with .= minimizes allocations by not creating small/temporary arrays
        @views for j in p.connections[i]
            du[2*i:2*i+1] .= du[2*i:2*i+1] .+ τ*spring_force(u[2*i:2*i+1], u[2*j:2*j+1], p.springs)
        end
    end
end


"""
    Water!(du, u, p::RaftParameters, t)

Compute the right-hand-side of the differential equation controlling the motion of a raft particles 
whose velocities are equal to `u = (1 - α)v_water + α v_wind`.

The parameters `p` are given by [`RaftParameters`](@ref), but only `p.α` and `p.gd_model` are used.

This function is equivalent to [`Raft!`](@ref) in the case where `τ = 0`.
"""
function Water!(du, u, p::RaftParameters, t)
    du[1] = p.gd_model.dSdt(u, t)

    α = p.clumps.α

    for i = 1:floor(Int64, length(u)/2)
        x, y = u[2*i:2*i+1]

        du[2*i] = u_x(x, y, t, α)
        du[2*i+1] = u_y(x, y, t, α)
    end
end
