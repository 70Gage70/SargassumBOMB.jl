using LinearAlgebra: ⋅

include("interpolant-constructors.jl")

########################################################

"""
    MaterialDerivativeX(vector_field, x, y, t)

Compute the material derivative of the x component of the vector field `vector_field` at coordinates `(x, y, t)`.
"""
function MaterialDerivativeX(vector_field::VectorField2DInterpolantEQR, x::Real, y::Real, t::Real)
    return gradient(vector_field.u, x, y, t) ⋅ [vector_field.u(x, y, t), vector_field.v(x, y, t), 1.0]
end

"""
    MaterialDerivativeY(vector_field, x, y, t)

Compute the material derivative of the y component of the vector field `vector_field` at coordinates `(x, y, t)`.
"""
function MaterialDerivativeY(vector_field::VectorField2DInterpolantEQR, x::Real, y::Real, t::Real)
    return gradient(vector_field.v, x, y, t) ⋅ [vector_field.u(x, y, t), vector_field.v(x, y, t), 1.0]
end

"""
    Vorticity(vector_field, x, y, t)

Compute the vorticity of the vector field `vector_field` at coordinates `(x, y, t)`.
"""
function Vorticity(vector_field::VectorField2DInterpolantEQR, x::Real, y::Real, t::Real)
    return gradient(vector_field.v, x, y, t)[1] - gradient(vector_field.u, x, y, t)[2]
end