"""
    struct ClumpParameters{T}

A container for the high-level parameters of the BOM equations.

### Fields
- `α` []: The fraction of the wind field acting on the particle.
- `τ` [d]: Measures the inertial response time of the medium to the particle
- `R` []: A geometric parameter.
- `f` [1/d]: The Coriolis parameter in the σ plane.
- `σ` []: The Stokes drift parameter; this applies an additional fraction of the Stokes drift to the water velocity 
    component of the particle.
"""
struct ClumpParameters{T<:Real}
    α::T
    τ::T
    R::T
    f::T
    σ::T
end

"""
    ClumpParameters(; constants...)

Compute the parameters required for the BOM equations from physical constants.

### Arguments

- `δ` []: The bouyancy of the particle. Default: `1.25`.
- `a` [km]: The radius of the particle. Default: `1.0e-4`.
- `ρ` [kg/km^3]: The density of the water. Default: `1027.0e9`.
- `ρa` [kg/km^3]: The density of the air. Default: `1.2e9`.
- `ν` [km^2/d]: The viscosity of the water. Default: `8.64e-8`.
- `νa` [km^2/d]: The viscosity of the air. Default: `1.296e-6`.
- `Ω` [rad/d]: The angular velocity of the Earth. Default: `2π`.
- σ []: The Stokes drift parameter. Default: `0.0`.
"""
function ClumpParameters(;
    δ::Real = 1.25,
    a::Real = 1.0e-4,
    ρ::Real = 1027.0e9,
    ρa::Real = 1.2e9,
    ν::Real = 8.64e-8,
    νa::Real = 1.296e-6,
    Ω::Real = 2*π,
    σ::Real = 0.0)

    μ = ν * ρ
    μa = νa * ρa
    γ = μa/μ

    ψ = (im*sqrt(1 - (2/δ - 1)^2) + 2/δ - 1)^(1/3)
    Φ = real(im*sqrt(3)/2 * (1/ψ - ψ) - 1/2 * (1/ψ + ψ) + 1)
    Ψ = 1/π * acos(1 - Φ) - 1/π * (1 - Φ) * sqrt(1 - (1 - Φ)^2)

    α = γ*Ψ/(1 - (1 - γ)*Ψ)
    τ = (1 - Φ/6)/(1 - (1 - γ)*Ψ) * (a^2 * ρ / (3*μ*δ^4))
    R = (1 - Φ/2)/(1 - Φ/6)

    ϑ0 = EQR.x.lat0
    f = 2*Ω*sin(ϑ0*π/180)

    return ClumpParameters(α, τ, R, f, σ)
end

"""
    mutable struct RaftParameters{T, U, F, C, G, L}

A container for the parameters defining a raft. Each clump and spring are identical.

### Structure 

`RaftParameters` acts as the parameter container for [`Raft!`](@ref). The solution vector `u` is a vector of length `2n_clumps ` 
of the form `[x1, y1, x2, y2, ...]` giving the `x, y` coordinates of each clump.

### Fields
- `ics`: An [`InitialConditions`](@ref).
- `clumps`: The [`ClumpParameters`](@ref) shared by each clump in the raft.
- `springs`: A subtybe of [`AbstractSpring`](@ref).
- `n_clumps_tot`: An `Integer` equal to the total number of clumps that have ever existed (i.e. it is at least the number of clumps that exist at any specific time.)
- `connections`: A subtybe of [`AbstractConnections`](@ref).
- `loc2label`: A `Dict` such that `loc2label[t]` is itself a `Dict` mapping vector indices to the absolute label of the clump in that location at
the `i`th time step. For example, `loc2label[t0][j] = j` since, at the initial time `t0`, the `j`th location contains the `j`th clump. If 
clump 1 dies at some later time `t`, then `loc2label[t][1] = 2`, `loc2label[t][2] = 3` since every clump is shifted by one to the left.
- `gd_model`: A subtype of [`AbstractGrowthDeathModel`](@ref). 
- `land`:: A subtype of [`AbstractLand`](@ref).

### Constructors 

Use `RaftParameters(; ics, clumps, springs, connections, gd_model, land)`.
The quantities `n_clumps_tot` and `loc2label` are computed automatically.
"""
mutable struct RaftParameters{T<:Real, U<:Integer, S<:AbstractSpring, C<:AbstractConnections, G<:AbstractGrowthDeathModel, L<:AbstractLand}
    ics::InitialConditions{T}
    clumps::ClumpParameters{T}
    springs::S
    n_clumps_tot::U
    connections::C
    loc2label::Dict{T, Dict{U, U}}
    gd_model::G
    land::L

    function RaftParameters(;
        ics::InitialConditions{T},
        clumps::ClumpParameters{T},
        springs::S,
        connections::C,
        gd_model::G,
        land::L) where {T<:Real, S<:AbstractSpring, C<:AbstractConnections, G<:AbstractGrowthDeathModel, L<:AbstractLand}

        n_c = n_clumps(ics.ics)
        loc2label = Dict(ics.tspan[1] => Dict(i => i for i = 1:n_c))
        form_connections!(connections, ics.ics)

        return new{T, Int64, S, C, G, L}(ics, clumps, springs, n_c, connections, loc2label, gd_model, land)
    end
end