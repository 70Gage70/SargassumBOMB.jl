using DifferentialEquations
using ModelingToolkit

include("parameters.jl")
include("vector-field-methods.jl")

##################################################

@variables t
ddt = Differential(t)

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

@register_symbolic spring_force_x(x1, x2, y1, y2, parameters)
@register_symbolic spring_force_y(x1, x2, y1, y2, parameters)

function Clump(
    xy0::Vector{<:Real},
    clump_parameters::BOMParameters;
    name::Symbol,
    forced::Bool = false)

    @assert length(xy0) == 2 "Must provide `xy0` as a vector with two elements."

    ps = @parameters α = clump_parameters.α τ = clump_parameters.τ R = clump_parameters.R f = clump_parameters.f
    @variables x(t) = xy0[1] y(t) = xy0[2]
    @variables Fx(t) Fy(t)

    eqs = [
    ddt(x) ~ u_x(x, y, t, α) + τ * (
        R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α) + Fx
    ),
    ddt(y) ~ u_y(x, y, t, α) + τ * (
        R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α) + Fy
    )
    ]

    if !forced
        forces = [Fx ~ 0, Fy ~ 0]
        eqs = [eqs ; forces]
    end

    return ODESystem(eqs, t, [x, y, Fx, Fy], ps; name)
end


function Net(
    xy0::Matrix{<:Real},
    clump_parameters::Vector{<:BOMParameters},
    spring_parameters::Matrix{<:SpringParameters};
    name::Symbol)

    @assert size(xy0, 2) == 2 "`xy0` must be an N x 2 matrix."
    @assert size(xy0, 1) == length(clump_parameters) == size(spring_parameters, 1) == size(spring_parameters, 2) "The dimensions of the input variables must be consistent."

    N_clumps = size(xy0, 1)
    @named clump 1:N_clumps i -> Clump(xy0[i,:], clump_parameters[i], forced = true)

    forces_x = [clump[i].Fx ~ sum([spring_force_x(clump[i].x, clump[j].x, clump[i].y, clump[j].y, spring_parameters[i, j])] for j = 1:N_clumps if j != i)[1]
        for i = 1:N_clumps
    ]

    forces_y = [clump[i].Fy ~ sum([spring_force_y(clump[i].x, clump[j].x, clump[i].y, clump[j].y, spring_parameters[i, j])] for j = 1:N_clumps if j != i)[1]
        for i = 1:N_clumps
    ]

    eqs = [forces_x ; forces_y]

    return compose(ODESystem(eqs, t; name = name), clump[1:N_clumps]...)
end

# xy01 = [1.0, 2.0]
# xy0 = [1.0 2.0 ; 3.0 4.0]
# ref = EquirectangularReference(lon0 = -75.0, lat0 = 10.0);
# cp1 = BOMParameters(ref);
# cps = [cp1, cp1]

# k_const(d) = 3
# sp = [SpringParameters(k_const, 0) SpringParameters(k_const, 4) ; SpringParameters(k_const, 6) SpringParameters(k_const, 0)];

# @named clump_no_force = Clump(xy01, cp1)
# @named clump_with_force = Clump(xy01, cp1, forced = true)
# @named net = Net(xy0, cps, sp)