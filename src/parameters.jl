include("coordinates.jl")

#################################

"""
    struct ClumpParameters{T}

A container for the high-level parameters of the BOM equations.

### Fields
- `ref`: The `EquirectangularReference` with which the projection is defined.
- `α` []: The fraction of the wind field acting on the particle.
- `τ` [d]: Measures the inertial response time of the medium to the particle
- `R` []: A geometric parameter.
- `f` [1/d]: The Coriolis parameter in the β plane.
"""
struct ClumpParameters{T<:Real}
    ref::EquirectangularReference{T}
    α::T
    τ::T
    R::T
    f::T
end

"""
    ClumpParameters(;constants...)

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
function ClumpParameters(
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

    return ClumpParameters(ref, α, τ, R, f)
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

"""
    spring_force(xy1, xy2, parameters)

Calculate the x and y components of the force on a point particle with coordinates `xy1` which is attached by a spring defined by `parameters` to another point particle with coordinates `xy2`.
"""
function spring_force(xy1::Vector{<:Real}, xy2::Vector{<:Real}, parameters::SpringParameters)
    d = norm(xy1 - xy2)
    if isapprox(d, 0.0)
        return 0.0
    else
        return parameters.k(d)*(parameters.L/norm(xy1 - xy2) - 1)*(xy1 - xy2)
    end
end

# """
#     spring_force(xy1, xy2, parameters)

# Calculate the x and y components of the force on a point particle with coordinates `xy1` which is attached by a spring defined by `parameters` to another point particle with coordinates `xy2`.
# """
# function spring_force(xy1::Vector{<:Real}, xy2::Vector{<:Real}, parameters::AbstractDict)
#     k, L = (parameters[:k], parameters[:L])
#     d = norm(xy1 - xy2)
#     return k(d)*(L/d - 1)*(xy1 - xy2)
# end

"""
    struct RaftParameters{C, S}

A container for the parameters defining a raft.

### Fields
- `xy0`: A `N x 2` matrix such that `xy0[i,:]` gives the `[x, y]` coordinates of the `i`th clump.
- `clumps`: A vector of [`ClumpParameters`](@ref) such that `clumps[i]` gives the parameters of the `i`th clump.
- `springs`: A matrix of [`SpringParameters`](@ref) such that `springs[i, j]` gives the parameters of the spring joining the `i`th and `j`th clump.
"""
struct RaftParameters{T<:Matrix{<:Real}, C<:Vector{<:ClumpParameters}, S<:Matrix{<:SpringParameters}}
    xy0::T
    clumps::C
    springs::S
end

"""
    RectangularRaft(x_range, y_range, clump_parameters, spring_parameters; network_type, name)

Construct [`RaftParameters`](@ref) in a rectangular arrangement such that each clump and spring have the same parameters.

### Arguments

- `x_range` [km]: A range which gives the x coordinates of the clumps in the raft. Should be increasing.
- `y_range` [km]: A range which gives the y coordinates of the clumps in the raft. Should be increasing.
- `clump_parameters`: The [`ClumpParameters`](@ref) shared by each clump.
- `spring_parameters`: The [`spring_parameters`](@ref) shared by each spring.

### Optional Arguments

- `network_type`: How the springs are conencted in the raft.
    - `"nearest"`: The default value. Each clump is connected to its perpendicular neighbors.
    - `"full"`: Each clump is connected to each other clump.
    - `"none"`: No clumps are connected.
"""
function RectangularRaftParameters(
    x_range::AbstractRange{<:Real}, 
    y_range::AbstractRange{<:Real},
    clump_parameters::ClumpParameters, 
    spring_parameters::SpringParameters; 
    network_type::String = "nearest")

    @assert network_type in ["nearest", "full", "none"] "`network_type` not recognized."
    @assert step(x_range) > 0 "x range should be increasing."
    @assert step(y_range) > 0 "y range should be increasing."

    # a rectangular mesh
    network = reverse.(collect(Iterators.product(reverse(y_range), x_range))) # reverse so that the first row has the largest y
    n_col = length(x_range)
    n_row = length(y_range)
    N_clumps = length(network)

    xy0 = Matrix{Float64}(undef, N_clumps, 2)
    clump_parameters_raft = Vector{ClumpParameters}(undef, N_clumps)
    spring_parameters_raft = Matrix{SpringParameters}(undef, N_clumps, N_clumps)

    # these arrays are all constructed in dictionary order (across rows, then down columns)
    n(i, j) = (i - 1) * n_col + j

    for i = 1:n_row
        for j = 1:n_col
            # linearized initial conditions
            xy0[n(i, j), :] .= network[i, j] 

            # all clumps have the same parameters
            clump_parameters_raft[n(i, j)] = clump_parameters

            # identify the appropriate connections
            if network_type == "nearest"
                connections = filter(idx -> (1 <= idx[1] <= n_row) && (1 <= idx[2] <= n_col), [(i-1, j), (i+1, j), (i, j-1), (i, j+1)])
            elseif network_type == "full"
                connections = [(a, b) for a = 1:n_row for b = 1:n_col if (a, b) != (i, j)]
            elseif network_type == "none"
                connections = []
            end
            
            connections = map(x->n(x...), connections)

            # all connected springs have the same parameters, `spring_parameters`
            # all unconnected springs have SpringParameters(k->0, 0) (effectively no spring)
            spring_parameters_raft[n(i, j), connections] .= spring_parameters
            spring_parameters_raft[n(i, j), setdiff(1:N_clumps, connections)] .= SpringParameters(k->0, 0)
        end
    end

    return RaftParameters(xy0, clump_parameters_raft, spring_parameters_raft)
end