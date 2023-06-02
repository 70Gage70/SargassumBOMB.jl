"""
    BOM_parameters(;constants...)

Compute the parameters required for the BOM equations from physical constants.

## Arguments [units]
- δ []: The bouancy of the particle.
- a [km]: The radius of the particle.
- ρ [kg/km^3]: The density of the water.
- ρa [kg/km^3]: The density of the air.
- ν [km^2/d]: The viscosity of the water.
- νa [km^2/d]: The viscosity of the air.
- Ω [rad/d]: The angular velocity of the Earth.
- ϑ0 [deg]: A reference latitude for the β plane.
"""
function BOM_parameters(;
    δ::Real = 1.25,
    a::Real = 1.0e-4,
    ρ::Real = 1027.0e9,
    ρa::Real = 1.2e9,
    ν::Real = 8.64e-8,
    νa::Real = 1.296e-6,
    Ω::Real = 2*π,
    ϑ0::Real = 10.0)

    μ = ν * ρ
    μa = νa * ρa
    γ = μa/μ

    ψ = (im*sqrt(1 - (2/δ - 1)^2) + 2/δ - 1)^(1/3)
    Φ = real(im*sqrt(3)/2 * (1/ψ - ψ) - 1/2 * (1/ψ + ψ) + 1)
    Ψ = 1/π * acos(1 - Φ) - 1/π * (1 - Φ) * sqrt(1 - (1 - Φ)^2)

    α = γ*Ψ/(1 - (1 - γ)*Ψ)
    τ = (1 - Φ/6)/(1 - (1 - γ)*Ψ) * (a^2 * ρ / (3*μ*δ^4))
    R = (1 - Φ/2)/(1 - Φ/6)

    f = 2*Ω*sin(ϑ0*π/180)

    return (α, τ, R, f)
end