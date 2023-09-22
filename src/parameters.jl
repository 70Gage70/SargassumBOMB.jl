using LinearAlgebra: norm
using SargassumFromAFAI
using NearestNeighbors

include("coordinates.jl")

#################################

"""
    abstract type AbstractGrowthDeathModel

The abstract type for growth and death models.
"""
abstract type AbstractGrowthDeathModel end 

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

- `ref`: The `EquirectangularReference` with which the projection is defined.
- `δ` []: The bouancy of the particle.
- `a` [km]: The radius of the particle.
- `ρ` [kg/km^3]: The density of the water.
- `ρa` [kg/km^3]: The density of the air.
- `ν` [km^2/d]: The viscosity of the water.
- `νa` [km^2/d]: The viscosity of the air.
- `Ω` [rad/d]: The angular velocity of the Earth.
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
- `k` [kg/d^2]: A scalar function of one variable which represents the stiffness of the spring. Recover a spring constant by providing, e.g. `k -> 5`.
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
    spring_force(xy1, xy2, parameters)

Calculate the x and y components of the force on a point particle with coordinates `xy1` which is attached by a spring defined by `parameters` to another point particle with coordinates `xy2`.
"""
function spring_force(xy1::Vector{<:Real}, xy2::Vector{<:Real}, parameters::SpringParameters)
    d = norm(xy1 - xy2)
    return parameters.k(d)*(parameters.L/d - 1)*(xy1 - xy2)
end

"""
    mutable struct RaftParameters{T, U, F}

A container for the parameters defining a raft. Each clump and spring are identical.

### Structure 

`RaftParameters` acts as the parameter container for [`Raft!`](@ref). The solution vector `u` is a vector of length `2n_clumps + 1` 
such that `u[1]` is an "amount" parameter which controls the growth and death of clumps by biophysical effects. Then, 
`u[2*i:2*i+1]` for `i = 1:n_clumps` gives the `[x, y]` coordinates of the clump in position `i`.

### Fields
- `ics`: A `Vector` such that `ics[1] = n_clumps` and `ics[2:2*n_clumps+1]` represents the `[x, y]` coordinates of each clump.
- `clumps`: The [`ClumpParameters`](@ref) shared by each clump in the raft.
- `springs`: The [`SpringParameters`](@ref) shared by each spring joining the clumps.
- `n_clumps_tot`: An `Integer` equal to the total number of clumps that have ever existed (i.e. it is at least the number of clumps that exist at any specific time.)
- `connections`: A `Dict` such that `connections[idx]` is a vector of indices `idx'` where a spring is connected between clumps `idx` and `idx'`. This should be updated in-place as clumps grow and die, i.e. `connections` only shows the current connections and refers to vector indices, not absolute clump labels.
- `loc2label`: A `Dict` such that `loc2label[t]` is itself a `Dict` mapping vector indices to the absolute label of the clump in that location at
the `i`th time step. For example, `loc2label[t0][j] = j` since, at the initial time `t0`, the `j`th location contains the `j`th clump. If 
clump 1 dies at some later time `t`, then `loc2label[t][1] = 2`, `loc2label[t][2] = 3` since every clump is shifted by one to the left.
- `gd_model`: A subtype of `AbstractGrowthDeathModel`. Must implement `growths`, `deaths` and `dSdt` callable at the solution vector `u`.
"""
mutable struct RaftParameters{T<:Real, U<:Integer, F<:Function, G<:AbstractGrowthDeathModel}
    ics::Vector{T}
    clumps::ClumpParameters{T}
    springs::SpringParameters{F, T}
    n_clumps_tot::U
    connections::Dict{U, Vector{U}}
    loc2label::Dict{T, Dict{U, U}}
    gd_model::G
end

function Base.show(io::IO, x::RaftParameters)
    print(io, "RaftParameters[")
    show(io, Integer(x.ics[1]))
    print(io, " Clumps]")
end

"""
    RectangularRaftParameters(x_range, y_range, clump_parameters, spring_parameters, t0, network_type, gd_model)

Construct [`RaftParameters`](@ref) in a rectangular arrangement.

Use [`OneClumpRaftParameters`](@ref) to construct a raft with single clump.

### Arguments

- `x_range`: A range which gives the x coordinates of the clumps in the raft. Should be increasing and in equirectangular coordinates, not longitudes.
- `y_range`: A range which gives the y coordinates of the clumps in the raft. Should be increasing and in equirectangular coordinates, not latitudes.
- `clump_parameters`: The [`ClumpParameters`](@ref) shared by each clump.
- `spring_k`: The spring function `k(x)` as in `F = - k(x)(x - L)`. The natural length of the spring `L` is set automatically according to the 
distances between adjacent clumps.
- `t0`: The initial time.
- `network_type`: How the springs are conencted in the raft.
    - `"nearest"`: The default value. Each clump is connected to its perpendicular neighbors.
    - `"full"`: Each clump is connected to each other clump.
    - `"none"`: No clumps are connected.
- `gd_model`: A subtype of `AbstractGrowthDeathModel`. Must implement `growths`, `deaths` and `dSdt` callable at the solution vector `u`.
"""
function RectangularRaftParameters(
    x_range::AbstractRange{<:Real}, 
    y_range::AbstractRange{<:Real},
    clump_parameters::ClumpParameters, 
    spring_k::Function,
    t0::Real, 
    network_type::String,
    gd_model::AbstractGrowthDeathModel)

    @assert network_type in ["nearest", "full", "none"] "`network_type` not recognized."
    @assert step(x_range) > 0 "x range should be inreasing."
    @assert step(y_range) > 0 "y range should be increasing."

    # a rectangular mesh
    network = reverse.(collect(Iterators.product(reverse(y_range), x_range))) # reverse so that the first row has the largest y
    network = reshape(stack(network, dims = 1), (size(network, 1), size(network, 2), 2)) # convert from matrix of tuples to array
    n_col = length(x_range)
    n_row = length(y_range)
    n_clumps = n_col *n_row
    connections = Dict{NTuple{2, Int64}, Vector{NTuple{2, Int64}}}()

    for i = 1:n_row, j = 1:n_col
        if network_type == "nearest"
            connections[(i, j)] = filter(idx -> (1 <= idx[1] <= n_row) && (1 <= idx[2] <= n_col), [(i-1, j), (i+1, j), (i, j-1), (i, j+1)])
        elseif network_type == "full"
            connections[(i, j)] = [(a, b) for a = 1:n_row for b = 1:n_col if (a, b) != (i, j)]
        elseif network_type == "none"
            connections[(i, j)] = Vector{NTuple{2, Int64}}()
        end
    end

    L = (step(x_range) + step(y_range))/2
    spring_parameters = SpringParameters(spring_k, L)

    # now flatten all quantities
    ics = [n_clumps]
    ics = vcat(ics, [network[i, j, k] for i = 1:n_row for j = 1:n_col for k = 1:2])

    n(i, j) = (i - 1) * n_col + j

    connections_flat = Dict{Int64, Vector{Int64}}()

    for key in keys(connections)
        connections_flat[n(key...)] = connections[key] .|> x -> n(x...)
    end

    loc2label = Dict(t0 => Dict(i => i for i = 1:n_clumps))

    return RaftParameters(ics, clump_parameters, spring_parameters, n_clumps, connections_flat, loc2label, gd_model)
end

"""
    OneClumpRaftParameters(x0, y0, clump_parameters, t0, gd_model)

Construct [`RaftParameters`](@ref) with a single clump.

### Arguments

- `x0`: The x coordinate of the clump.
- `y0`: The x coordinate of the clump.
- `clump_parameters`: The [`ClumpParameters`](@ref) of the clump.
- `t0`: The initial time.
- `gd_model`: A subtype of `AbstractGrowthDeathModel`. Must implement `growths`, `deaths` and `dSdt` callable at the solution vector `u`.
"""
function OneClumpRaftParameters(
    x0::Real, 
    y0::Real, 
    clump_parameters::ClumpParameters, 
    t0::Real, 
    gd_model::AbstractGrowthDeathModel)

    n_clumps = 1
    ics = [1.0, x0, y0]
    spring_parameters = SpringParameters(k -> 0.0, 1.0)
    connections = Dict(1 => Int64[])
    loc2label = Dict(t0 => Dict(1 => 1))

    return RaftParameters(ics, clump_parameters, spring_parameters, n_clumps, connections, loc2label, gd_model)
end

# dists = SargassumDistribution("/Users/gagebonner/Desktop/Repositories/SargassumFromAFAI.jl/data/dist-2018.nc")

"""
    DistributionRaftParameters(dist, number, sample_type, clump_parameters, spring_k, t0, network_type, gd_model)

Construct [`RaftParameters`](@ref) from a `SargassumDistribution`.

### Arguments 
- `dist`: A `SargassumDistribution`.
- `number`: The number of clumps to initialize; interactive with `sample_type`.
- `sample_type`: A `String` identifying one of the methods of assigning clump locations based on the distribution.
    - `"sample"`: A number `number` of samples are drawn from `dist`. Each sample is placed uniformly at random inside corresponding box.
    - `"sorted"`: Boxed are filled with one clump placed uniformly at random inside them, starting from the box with the highest concentration. If `number` 
                is greater than the total number of boxes, repeat the loop starting again from the highest concentration box.
    - `"uniform"`: Exactly one clump is placed in the center of each box with nonzero concentration. Note that this ignores `number`.
- `clump_parameters`: The [`ClumpParameters`](@ref) shared by each clump.
- `spring_k`: The spring function `k(x)` as in `F = - k(x)(x - L)`. The natural length of the spring `L` is set automatically according to the 
distances between adjacent clumps.
- `t0`: The initial time.
- `network_type`: How the springs are conencted in the raft.
    - `"nearest"`: Each clump is connected to all neighbors within a radius of approximately two gridpoints.
    - `"full"`: Each clump is connected to each other clump.
    - `"none"`: No clumps are connected.
- `gd_model`: A subtype of `AbstractGrowthDeathModel`. Must implement `growths`, `deaths` and `dSdt` callable at the solution vector `u`.
"""
function DistributionRaftParameters(
    dist::SargassumDistribution, 
    number::Integer,
    sample_type::String,
    clump_parameters::ClumpParameters, 
    spring_k::Function,
    t0::Real, 
    network_type::String,
    gd_model::AbstractGrowthDeathModel)

    @assert sample_type in ["sample", "sorted", "uniform"] "`sample_type` not recognized."
    @assert network_type in ["nearest", "full", "none"] "`network_type` not recognized."

    sarg = dist.sargassum
    lons = dist.lon
    lats = dist.lat

    δ_x = abs(lons[2] - lons[1])
    δ_y = abs(lats[2] - lats[1])

    xy0 = Float64[]

    if sample_type == "sample"
        pts = Iterators.product(lons, lats) |> x -> reduce(vcat, collect(x))
        wts = sarg |> x -> reduce(vcat, x)
        samples = sample(pts, Weights(wts), number)

        for samp in samples
            push!(xy0, rand(Uniform(samp[1] - δ_x/2, samp[1] + δ_x/2)))
            push!(xy0, rand(Uniform(samp[2] - δ_y/2, samp[2] + δ_y/2)))
        end
    elseif sample_type == "sorted"
        pts = Iterators.product(lons, lats) |> x -> reduce(vcat, collect(x))
        wts = sarg |> x -> reduce(vcat, x)

        idx = findall(x -> x > 0, wts)
        pts = pts[idx]
        wts = wts[idx]

        pts = pts[sortperm(wts, rev = true)]

        n_c = 1
        while n_c <= number
            idx = n_c == length(pts) ? length(pts) : mod(n_c, length(pts))
            pt = pts[idx]

            push!(xy0, rand(Uniform(pt[1] - δ_x/2, pt[1] + δ_x/2)))
            push!(xy0, rand(Uniform(pt[2] - δ_y/2, pt[2] + δ_y/2)))       

            n_c = n_c + 1
        end
    elseif sample_type == "uniform"
        pts = Iterators.product(lons, lats) |> x -> reduce(vcat, collect(x))
        wts = sarg |> x -> reduce(vcat, x)

        idx = findall(x -> x > 0, wts)
        pts = pts[idx]

        for pt in pts
            push!(xy0, pt[1])
            push!(xy0, pt[2])  
        end
    end

    n_clumps = length(xy0)/2
    pushfirst!(xy0, n_clumps)
    n_clumps = Int64(n_clumps)

    L = (δ_x + δ_y)/2
    spring_parameters = SpringParameters(spring_k, L)

    if network_type == "full"
        connections = Dict(1:n_clumps .=> [[j for j = 1:n_clumps if j != i] for i = 1:n_clumps])
    elseif network_type == "none"
        connections = Dict(1:n_clumps .=> [Int64[] for i = 1:n_clumps])
    elseif network_type == "nearest"
        data = reshape(xy0[2:end], 2, n_clumps)
        balltree = BallTree(data)
        idx = inrange(balltree, data, 2*L, true)
        idx = [filter(x -> x != i, idx[i]) for i = 1:length(idx)]
        connections = Dict(1:n_clumps .=> idx)
    end

    loc2label = Dict(t0 => Dict(i => i for i = 1:n_clumps))

    return RaftParameters(xy0, clump_parameters, spring_parameters, n_clumps, connections, loc2label, gd_model)
end

