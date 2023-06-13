include("coordinates.jl")

#################################

"""
    struct BOMParameters{T}

A container for the high-level parameters of the BOM equations.

### Fields
- `ref`: The `EquirectangularReference` with which the projection is defined.
- `α` []: The fraction of the wind field acting on the particle.
- `τ` [d]: Measures the inertial response time of the medium to the particle
- `R` []: A geometric parameter.
- `f` [1/d]: The Coriolis parameter in the β plane.
"""
struct BOMParameters{T<:Real}
    ref::EquirectangularReference{T}
    α::T
    τ::T
    R::T
    f::T
end

"""
    BOMParameters(;constants...)

Compute the parameters required for the BOM equations from physical constants.

### Arguments [units]

- `δ` []: The bouancy of the particle.
- `a` [km]: The radius of the particle.
- `ρ` [kg/km^3]: The density of the water.
- `ρa` [kg/km^3]: The density of the air.
- `ν` [km^2/d]: The viscosity of the water.
- `νa` [km^2/d]: The viscosity of the air.
- `Ω` [rad/d]: The angular velocity of the Earth.
- `ref`: The `EquirectangularReference` with which the projection is defined.
"""
function BOMParameters(
    ref::EquirectangularReference;
    δ::Real = 1.25,
    a::Real = 1.0e-4,
    ρ::Real = 1027.0e9,
    ρa::Real = 1.2e9,
    ν::Real = 8.64e-8,
    νa::Real = 1.296e-6,
    Ω::Real = 2*π)

    μ = ν * ρ
    μa = νa * ρa
    γ = μa/μ

    ψ = (im*sqrt(1 - (2/δ - 1)^2) + 2/δ - 1)^(1/3)
    Φ = real(im*sqrt(3)/2 * (1/ψ - ψ) - 1/2 * (1/ψ + ψ) + 1)
    Ψ = 1/π * acos(1 - Φ) - 1/π * (1 - Φ) * sqrt(1 - (1 - Φ)^2)

    α = γ*Ψ/(1 - (1 - γ)*Ψ)
    τ = (1 - Φ/6)/(1 - (1 - γ)*Ψ) * (a^2 * ρ / (3*μ*δ^4))
    R = (1 - Φ/2)/(1 - Φ/6)

    ϑ0 = ref.lat0
    f = 2*Ω*sin(ϑ0*π/180)

    return BOMParameters(ref, α, τ, R, f)
end

"""
    SpringParameters{T}

A container for the parameters defining a spring.
   
### Fields
- `k` [kg/d^2]: A scalar function of one variable which represents the stiffness of the spring. Recover a spring constant by providing, e.g. k(d) =  5.
- `L` [km]: The natural length of the spring.
"""
struct SpringParameters{F<:Function, T<:Real}
    k::F
    L::T
end

"""
    spring_force_x(x1, x2, y1, y2, parameters)

Calculate the x component of the force on a point particle with coordinates `(x1, y1)` which is attached by a spring defined by `parameters` to another point particle with coordinates `(x2, y2)`.
"""
function spring_force_x(x1::Real, x2::Real, y1::Real, y2::Real, parameters::SpringParameters)
    k, L = (parameters.k, parameters.L)
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    return k(d)*(x1 - x2)*(L/d - 1)
end

"""
    spring_force_y(x1, x2, y1, y2, parameters)

Calculate the y component of the force on a point particle with coordinates `(x1, y1)` which is attached by a spring defined by `parameters` to another point particle with coordinates `(x2, y2)`.
"""
function spring_force_y(x1::Real, x2::Real, y1::Real, y2::Real, parameters::SpringParameters)
    k, L = (parameters.k, parameters.L)
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    return k(d)*(y1 - y2)*(L/d - 1)
end