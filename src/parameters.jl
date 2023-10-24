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
- `β` []: The Stokes drift parameter; this applies an additional fraction of the wind to the particle.
"""
struct ClumpParameters{T<:Real}
    ref::EquirectangularReference{T}
    α::T
    τ::T
    R::T
    f::T
    β::T
end

"""
    ClumpParameters(ref; constants...)

Compute the parameters required for the BOM equations from physical constants.

### Arguments [units]

- `ref`: The `EquirectangularReference` with which the projection is defined.
- `δ` []: The bouancy of the particle. Default: `1.25`.
- `a` [km]: The radius of the particle. Default: `1.0e-4`.
- `ρ` [kg/km^3]: The density of the water. Default: `1027.0e9`.
- `ρa` [kg/km^3]: The density of the air. Default: `1.2e9`.
- `ν` [km^2/d]: The viscosity of the water. Default: `8.64e-8`.
- `νa` [km^2/d]: The viscosity of the air. Default: `1.296e-6`.
- `Ω` [rad/d]: The angular velocity of the Earth. Default: `2π`.
- β []: The Stokes drift parameter. Default: `0.0`.
"""
function ClumpParameters(
    ref::EquirectangularReference;
    δ::Real = 1.25,
    a::Real = 1.0e-4,
    ρ::Real = 1027.0e9,
    ρa::Real = 1.2e9,
    ν::Real = 8.64e-8,
    νa::Real = 1.296e-6,
    Ω::Real = 2*π,
    β::Real = 0.0)

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

    return ClumpParameters(ref, α, τ, R, f, β)
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

function spring_force(
    xy1::SubArray{T, 1, Vector{T}, Tuple{UnitRange{R}}},
    xy2::SubArray{T, 1, Vector{T}, Tuple{UnitRange{R}}},
    parameters::SpringParameters) where {T <:Real, R<:Integer}
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
- `gd_model`: A subtype of `AbstractGrowthDeathModel`. Must implement `growths`, `deaths` and `dSdt` callable at the solution vector and current time, 
i.e. it must be a function `dSdt(u, t)`.

### Constructors 

Use `RaftParameters(; ics, clumps, springs, connections, t0, gd_model)`. The arguments `clumps`, `springs` and 
`gd_model` are passed directly to the struct.

The [`initial_conditions`](@ref) methods help with the construction of `ics` and the function [`form_connections`](@ref) 
helps with the construction of `connections`.
"""
mutable struct RaftParameters{T<:Real, U<:Integer, F<:Function, G<:AbstractGrowthDeathModel}
    ics::Vector{T}
    clumps::ClumpParameters{T}
    springs::SpringParameters{F, T}
    n_clumps_tot::U
    connections::Dict{U, Vector{U}}
    loc2label::Dict{T, Dict{U, U}}
    gd_model::G

    function RaftParameters(;
        ics::Vector{T},
        clumps::ClumpParameters{T},
        springs::SpringParameters{F, T},
        connections::Dict{U, Vector{U}},
        t0::T,
        gd_model::G) where {T<:Real, U<:Integer, F<:Function, G<:AbstractGrowthDeathModel}

        n_clumps = Int64((length(ics) - 1)/2)
        loc2label = Dict(t0 => Dict(i => i for i = 1:n_clumps))

        return new{T, U, F, G}(ics, clumps, springs, n_clumps, connections, loc2label, gd_model)
    end
end

function Base.show(io::IO, x::RaftParameters)
    print(io, "RaftParameters[")
    show(io, Integer(x.ics[1]))
    print(io, " Clumps]")
end

"""
    initial_conditions(xy0; ref = nothing)

Construct initial conditions suitable for use in `RaftParameters.ics` from a list of coordinates of the form 
`[x1, y1, x2, y2 ..., xN, yN]`. These should be equirectangular coordinates; if `ref` is provided, the coordinates 
are converted from spherical coordinates.

Can be applied as `initial_conditions(x_range, y_range; ref = nothing)` to generate clumps in a rectangular arrangement.

Can be applied as `initial_conditions(x0, y0; ref = nothing)` for a single clump with coordinates `(x0, y0)`.
"""
function initial_conditions(xy0::Vector{<:Real}; ref::Union{Nothing, EquirectangularReference} = nothing)
    if ref !== nothing
        ics = sph2xy(xy0, ref)
    else
        ics = deepcopy(xy0)
    end

    pushfirst!(ics, length(xy0)/2)
    return ics
end

function initial_conditions(
    x_range::AbstractRange{T}, 
    y_range::AbstractRange{T}; 
    ref::Union{Nothing, EquirectangularReference} = nothing) where {T<:Real}

    if ref !== nothing
        ics_x, ics_y = sph2xy(x_range, y_range, ref)
    else
        ics_x, ics_y = x_range, y_range
    end

    ics = T[length(x_range)*length(y_range)]

    for x in ics_x, y in ics_y
        push!(ics, x, y)
    end

    return ics
end

function initial_conditions(x0::Real, y0::Real; ref::Union{Nothing, EquirectangularReference} = nothing)
    if ref !== nothing
        ics = sph2xy(x0, y0, ref)
    else
        ics = [x0, y0]
    end

    pushfirst!(ics, 1)
    return ics
end

"""
    initial_conditions(dist, number, sample_type, ref)

Construct [`RaftParameters`](@ref) from a `SargassumDistribution`.

### Arguments 
- `dist`: A `SargassumDistribution`.
- `weeks`: A `Vector{<:Integer}` giving the weeks of the month to consider. Each entry should be between 1 and 4 and appear only once.
- `number`: The number of clumps to initialize; interactive with `sample_type`.
- `sample_type`: A `String` identifying one of the methods of assigning clump locations based on the distribution.
    - `"sample"`: A number `number` of samples are drawn from `dist`. Each sample is placed uniformly at random inside the corresponding box.
    - `"sorted"`: Boxed are filled with one clump placed uniformly at random inside them, starting from the box with the highest concentration. If `number` 
                is greater than the total number of boxes, repeat the loop starting again from the highest concentration box.
    - `"uniform"`: Exactly one clump is placed in the center of each box with nonzero concentration. Note that this ignores `number`.
- `ref`: An [`EquirectangularReference`](@ref). A `SargassumDistribution` has fields `lon` and `lat`, so this is necessary to 
                covert these to equirectangular coordinates.
"""
function initial_conditions(
    dist::SargassumDistribution, 
    weeks::Vector{<:Integer},
    number::Integer,
    sample_type::String,
    ref::EquirectangularReference)

    @assert sample_type in ["sample", "sorted", "uniform"] "`sample_type` not recognized."
    @assert length(weeks) > 0 "At least one week must be selected."
    @assert all(map(x -> 1 <= x <= 4, weeks)) "Each entry of `weeks` must be between 1 and 4."
    @assert allunique(weeks) "Each week should only appear once."

    sarg = sum(dist.sargassum[:,:,week] for week in weeks)
    lons = dist.lon
    lats = dist.lat

    δ_x = abs(lons[2] - lons[1])
    δ_y = abs(lats[2] - lats[1])

    xy0 = Float64[]

    if sample_type == "sample"
        pts = Iterators.product(lons, lats) |> x -> reduce(vcat, collect(x))
        wts = sarg |> x -> reduce(vcat, x)
        samples = Distributions.sample(pts, Weights(wts), number)

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
            idx = mod(n_c, length(pts)) == 0 ? length(pts) : mod(n_c, length(pts))
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

    xy0 = sph2xy(xy0, ref)
    n_clumps = length(xy0)/2
    pushfirst!(xy0, n_clumps)
    
    return xy0
end

"""
    form_connections(ics, network_type; neighbor_parameter = nothing)

Construct connections between clumps suitable for use in `RaftParameters.connections`.

### Arguments 

- `ics`: A `Vector` of initial conditions as defined by [`RaftParameters`](@ref). Note that a solution vector of a [`Raft!`](@ref) problem 
    can actually be provided at any time to form connections dynamically as with [`cb_connections`](@ref).
- `network_type`: A `String` identifying how the connections should be consructed.
    - `"full"`: Every clump is connected to every other clump.
    - `"none"`: No clumps are connected.
    - `"radius"`: Each clump is connected to clumps within a distance `neighbor_parameter`. If `neighbor_parameter` is not provided, 
                    the distance is taken to be the mean value of all pairwise clump distances in `ics`.
    - `"nearest"`: Each clump is connected to its `neighbor_parameter` nearest neighbors (not including itself). If `neighbor_parameter` 
                    is not provided, each clump is connected to its nearest neighbor, i.e. `neighbor_parameter = 1`.

### Optional Arguments

- `neighbor_parameter`: A parameter controlling the connections for the `"radius"` and `"nearest"` options of `network_type`; see above.
"""
function form_connections(
    ics::Vector{<:Real}, 
    network_type::String; 
    neighbor_parameter::Union{Nothing, Real} = nothing)
    @assert network_type in ["full", "none", "radius", "nearest"] "`network_type` not recognized."

    n_clumps = Int64((length(ics) - 1)/2)

    if network_type == "full"
        connections = Dict(1:n_clumps .=> [[j for j = 1:n_clumps if j != i] for i = 1:n_clumps])
    elseif network_type == "none"
        connections = Dict(1:n_clumps .=> [Int64[] for i = 1:n_clumps])
    elseif network_type == "radius"
        data = reshape(@view(ics[2:end]), 2, n_clumps)

        if neighbor_parameter === nothing
            dists = Float64[]
            for i = 1:size(data, 2) - 1
                for j = i+1:size(data, 2)
                    push!(dists, norm(data[:,i] - data[:, j]))
                end
            end

            radius = mean(dists)
        else
            radius = neighbor_parameter
        end
        
        balltree = BallTree(data)
        idx = inrange(balltree, data, radius, true)
        idx = [filter(x -> x != i, idx[i]) for i = 1:length(idx)]
        connections = Dict(1:n_clumps .=> idx)
    elseif network_type == "nearest"
        # each point is considered one of its own nearest neighbors, so have to add 1
        if neighbor_parameter === nothing
            k = 2
        else
            k = neighbor_parameter + 1
        end
        
        k = min(n_clumps, k)

        data = reshape(@view(ics[2:end]), 2, n_clumps)
        kdtree = KDTree(data)
        idx, dists = knn(kdtree, data, k, true)
        idx = [filter(x -> x != i, idx[i]) for i = 1:length(idx)]
        connections = Dict(1:n_clumps .=> idx)
    end

    return connections
end