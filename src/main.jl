using DifferentialEquations
using ModelingToolkit

include("parameters.jl")
include("vector-fields.jl")

#############################################

@parameters α τ R f
@variables t x(t) y(t)

vfs = vector_fields()
v_x(x, y, t) = vfs[1](x, y, t)
v_y(x, y, t) = vfs[2](x, y, t)
Dv_xDt(x, y, t) = vfs[3](x, y, t)
Dv_yDt(x, y, t) = vfs[4](x, y, t)
u_x(x, y, t, α) = vfs[5](x, y, t, α)
u_y(x, y, t, α) = vfs[6](x, y, t, α)
Du_xDt(x, y, t, α) = vfs[7](x, y, t, α)
Du_yDt(x, y, t, α) = vfs[8](x, y, t, α)
ω(x, y, t) = vfs[9](x, y, t)

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