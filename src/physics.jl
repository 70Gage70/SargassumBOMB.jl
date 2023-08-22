using OrdinaryDiffEq
using JLD2

include(joinpath(@__DIR__, "parameters.jl"))
include(joinpath(@__DIR__, "..", "interpolants", "interpolant-constructors.jl"))
include(joinpath(@__DIR__, "..", "interpolants", "interpolant-derivatives.jl"))

########################################################################

# loading interpolants
itp_path = joinpath(@__DIR__, "..", "interpolants", "ocean-atmos")
isdefined(@__MODULE__, :water_itp) || (const water_itp = load(joinpath(itp_path, "water_itp.jld2"), "water_itp"))
isdefined(@__MODULE__, :wind_itp) || (const wind_itp = load(joinpath(itp_path, "wind_itp.jld2"), "wind_itp"))

# All the functions depending on wind and water vector fields.
v_x(x, y, t) = water_itp.u(x, y, t)
v_y(x, y, t) =  water_itp.v(x, y, t)
Dv_xDt(x, y, t) = MaterialDerivativeX(water_itp, x, y, t)
Dv_yDt(x, y, t) = MaterialDerivativeY(water_itp, x, y, t)
u_x(x, y, t, α) = (1 - α) * water_itp.u(x, y, t) + α * wind_itp.u(x, y, t)
u_y(x, y, t, α) = (1 - α) * water_itp.v(x, y, t) + α * wind_itp.v(x, y, t)
Du_xDt(x, y, t, α) = (1 - α) * MaterialDerivativeX(water_itp, x, y, t) + α * MaterialDerivativeX(wind_itp, x, y, t) 
Du_yDt(x, y, t, α) = (1 - α) * MaterialDerivativeY(water_itp, x, y, t) + α * MaterialDerivativeY(wind_itp, x, y, t) 
ω(x, y, t) = Vorticity(water_itp, x, y, t)

########################################################################
########################################################################
########################################################################

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


# have u[1:2] = (x, y)_1, u[3:4] = (x, y)_2 ... u[2*i-1:2*i] = (x, y)_i
function Raft!(du, u, p::RaftParameters, t)
    α, τ, R, f = p.clumps.α, p.clumps.τ, p.clumps.R, p.clumps.f

    for i = 1:Integer(length(u)/2)
        x, y = u[2*i-1:2*i]

        du[2*i-1] = u_x(x, y, t, α) + τ * (
            R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α)
        )
        du[2*i] = u_y(x, y, t, α) + τ * (
            R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α)
        )

        du[2*i-1:2*i] += τ*sum(spring_force(u[2*i-1:2*i], u[2*j-1:2*j], p.springs) for j in p.connections[i]; init = [0.0, 0.0])
    end
end

