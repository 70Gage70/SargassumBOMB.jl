v_x(x, y, t, σ) = WATER_ITP.x.fields[:u](x, y, t) + σ * STOKES_ITP.x.fields[:u](x, y, t)
v_y(x, y, t, σ) =  WATER_ITP.x.fields[:v](x, y, t) + σ * STOKES_ITP.x.fields[:v](x, y, t)
Dv_xDt(x, y, t, σ) = WATER_ITP.x.fields[:DDt_x](x, y, t) + σ * STOKES_ITP.x.fields[:DDt_x](x, y, t)
Dv_yDt(x, y, t, σ) = WATER_ITP.x.fields[:DDt_y](x, y, t) + σ * STOKES_ITP.x.fields[:DDt_y](x, y, t)
u_x(x, y, t, α, σ) = (1 - α) * v_x(x, y, t, σ) + α * WIND_ITP.x.fields[:u](x, y, t)
u_y(x, y, t, α, σ) = (1 - α) * v_y(x, y, t, σ) + α * WIND_ITP.x.fields[:v](x, y, t)
Du_xDt(x, y, t, α, σ) = (1 - α) * Dv_xDt(x, y, t, σ) + α * WIND_ITP.x.fields[:DDt_x](x, y, t)
Du_yDt(x, y, t, α, σ) = (1 - α) * Dv_yDt(x, y, t, σ) + α * WIND_ITP.x.fields[:DDt_y](x, y, t)
ω(x, y, t) = WATER_ITP.x.fields[:vorticity](x, y, t)

leeway_x(x, y, t, α, σ) = v_x(x, y, t, σ) + α * WIND_ITP.x.fields[:u](x, y, t)
leeway_y(x, y, t, α, σ) = v_y(x, y, t, σ) + α * WIND_ITP.x.fields[:v](x, y, t)
########################################################################

"""
    Raft!(du, u, p::RaftParameters, t)

Compute the right-hand-side of the differential equation controlling the motion of a raft with parameters 
given by [`RaftParameters`](@ref).

The solution vector `u` is a vector of length `2n_clumps + 1` such that `u[1]` is an "amount" parameter which 
controls the growth and death of clumps by biophysical effects. Then, `u[2*i:2*i+1]` for `i = 1:n_clumps` gives 
the `[x, y]` coordinates of the clump in position `i`.

For integrating using a leeway velocity, [`Leeway!`](@ref) should be used.
"""
function Raft!(du, u, p::RaftParameters, t)
    du[1] = p.gd_model.dSdt(u, t)

    α, τ, R, f, σ = p.clumps.α, p.clumps.τ, p.clumps.R, p.clumps.f, p.clumps.σ

    for i = 1:floor(Int64, length(u)/2)
        # note that x, y = u[2*i], u[2*i+1]

        du[2*i] = 
            u_x(u[2*i], u[2*i+1], t, α, σ) 
            + τ * (
            R*Dv_xDt(u[2*i], u[2*i+1], t, σ) 
            - R*(f + ω(u[2*i], u[2*i+1], t)/3)*v_y(u[2*i], u[2*i+1], t, σ) 
            - Du_xDt(u[2*i], u[2*i+1], t, α, σ) 
            + (f + R*ω(u[2*i], u[2*i+1], t)/3)*u_y(u[2*i], u[2*i+1], t, α, σ)
        )
        du[2*i+1] = 
            u_y(u[2*i], u[2*i+1], t, α, σ) 
            + τ * (
            R*Dv_yDt(u[2*i], u[2*i+1], t, σ) 
            + R*(f + ω(u[2*i], u[2*i+1], t)/3)*v_x(u[2*i], u[2*i+1], t, σ) 
            - Du_yDt(u[2*i], u[2*i+1], t, α, σ) 
            - (f + R*ω(u[2*i], u[2*i+1], t)/3)*u_x(u[2*i], u[2*i+1], t, α, σ)
        )

        for j in p.connections.connections[i]
            d = sqrt((u[2*i] - u[2*j])^2 + (u[2*i+1] - u[2*j+1])^2)
            fac = τ*(p.springs.k(d))*(p.springs.L/d - 1)

            du[2*i] += fac*(u[2*i] - u[2*j])
            du[2*i+1] += fac*(u[2*i+1] - u[2*j+1])
        end
    end
end


"""
    Leeway!(du, u, p::RaftParameters, t)

Compute the right-hand-side of the differential equation controlling the motion of raft particles 
whose velocities are equal to `u = v_water + σ v_stokes + α v_wind`.

The parameters `p` are given by [`RaftParameters`](@ref), but only `p.α, p.σ` and `p.gd_model` are used.
"""
function Leeway!(du, u, p::RaftParameters, t)
    du[1] = p.gd_model.dSdt(u, t)

    α, σ = p.clumps.α, p.clumps.σ

    for i = 1:floor(Int64, length(u)/2)
        x, y = u[2*i:2*i+1]

        du[2*i] = leeway_x(x, y, t, α, σ)
        du[2*i+1] = leeway_y(x, y, t, α, σ)
    end
end
