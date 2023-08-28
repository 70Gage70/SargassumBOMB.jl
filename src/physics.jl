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
v_x(x, y, t) = water_itp.fields[:u](x, y, t)
v_y(x, y, t) =  water_itp.fields[:v](x, y, t)
Dv_xDt(x, y, t) = MaterialDerivativeX(water_itp, x, y, t)
Dv_yDt(x, y, t) = MaterialDerivativeY(water_itp, x, y, t)
u_x(x, y, t, α) = (1 - α) * water_itp.fields[:u](x, y, t) + α * wind_itp.fields[:u](x, y, t)
u_y(x, y, t, α) = (1 - α) * water_itp.fields[:v](x, y, t) + α * wind_itp.fields[:v](x, y, t)
Du_xDt(x, y, t, α) = (1 - α) * MaterialDerivativeX(water_itp, x, y, t) + α * MaterialDerivativeX(wind_itp, x, y, t) 
Du_yDt(x, y, t, α) = (1 - α) * MaterialDerivativeY(water_itp, x, y, t) + α * MaterialDerivativeY(wind_itp, x, y, t) 
ω(x, y, t) = Vorticity(water_itp, x, y, t)

########################################################################
########################################################################
########################################################################

"""
    `Clump!`(du, u, p::ClumpParameters, t)

Compute the right-hand-side of the differential equation controlling the motion of a single clump with 
parameters given by [`ClumpParameters`](@ref).

The solution vector `u` is a 2d vector such that `u[1:2] = [x, y]`.

### Example 

```julia
x0, y0 = -64, 14
tspan = (0.0, 200.0)
xy0 = sph2xy(x0, y0, ref_itp) 

cp = ClumpParameters(ref_itp)
clump_prob = ODEProblem(Clump!, xy0, tspan, cp)
sol_clump = solve(clump_prob, Tsit5())
```
"""
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


"""
    `Raft!`(du, u, p::RaftParameters, t)

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
        x, y = u[2*i:2*i+1]

        du[2*i] = u_x(x, y, t, α) + τ * (
            R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α)
        )
        du[2*i+1] = u_y(x, y, t, α) + τ * (
            R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α)
        )

        du[2*i:2*i+1] += τ*sum(spring_force(u[2*i:2*i+1], u[2*j:2*j+1], p.springs) for j in p.connections[i]; init = [0.0, 0.0])
    end
end

