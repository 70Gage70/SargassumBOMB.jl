struct BOMParameters{T<:Real}
    α::T
    τ::T
    R::T
    f::T
end

function BOMParameters(;
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

    return BOMParameters(α, τ, R, f)
end

struct SpringParameters{F<:Function, T<:Real}
    k::F
    L::T
end

function Base.length(::SpringParameters)
    return 1
end

function Base.iterate(sp::SpringParameters)
    return (sp, nothing)
end

function Base.iterate(::SpringParameters, ::Nothing)
    return nothing
end

function Base.show(io::IO, x::SpringParameters)
    print(io, "SpringParameters[1->")
    show(io, x.k(1))
    print(io, ", ")
    show(io, length(x.L))
    print(io, "]")
end

function spring_force_x(x1::Real, x2::Real, y1::Real, y2::Real, parameters::SpringParameters)
    k, L = (parameters.k, parameters.L)
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    return k(d)*(x1 - x2)*(L/d - 1)
end

function spring_force_y(x1::Real, x2::Real, y1::Real, y2::Real, parameters::SpringParameters)
    k, L = (parameters.k, parameters.L)
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    return k(d)*(y1 - y2)*(L/d - 1)
end