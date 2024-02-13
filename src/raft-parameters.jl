"""
    struct ClumpParameters{T}

A container for the high-level parameters of the BOM equations.

### Fields
- `ref`: The `EquirectangularReference` with which the projection is defined.
- `α` []: The fraction of the wind field acting on the particle.
- `τ` [d]: Measures the inertial response time of the medium to the particle
- `R` []: A geometric parameter.
- `f` [1/d]: The Coriolis parameter in the σ plane.
- `σ` []: The Stokes drift parameter; this applies an additional fraction of the Stokes drift to the water velocity 
    component of the particle.
"""
struct ClumpParameters{T<:Real}
    ref::EquirectangularReference{T}
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

- `ref`: The `EquirectangularReference` with which the projection is defined. Default [`EQR_DEFAULT`](@ref).
- `δ` []: The bouancy of the particle. Default: `1.25`.
- `a` [km]: The radius of the particle. Default: `1.0e-4`.
- `ρ` [kg/km^3]: The density of the water. Default: `1027.0e9`.
- `ρa` [kg/km^3]: The density of the air. Default: `1.2e9`.
- `ν` [km^2/d]: The viscosity of the water. Default: `8.64e-8`.
- `νa` [km^2/d]: The viscosity of the air. Default: `1.296e-6`.
- `Ω` [rad/d]: The angular velocity of the Earth. Default: `2π`.
- σ []: The Stokes drift parameter. Default: `0.0`.
"""
function ClumpParameters(;
    ref::EquirectangularReference = EQR_DEFAULT,
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

    ϑ0 = ref.lat0
    f = 2*Ω*sin(ϑ0*π/180)

    return ClumpParameters(ref, α, τ, R, f, σ)
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
    show(io, x.L)
    print(io, "]")
end

"""
    ΔL(x_range, y_range; ref)

Compute a spring length from a rectangular arrangement of clumps provided by `x_range` and `y_range`. This is the distance between the centers of 
diagonally-adjacent gridpoints.

If `ref` is provided, the ranges are converted from spherical to equirectangular coordinates. Default `nothing`.
"""
function ΔL(x_range::AbstractRange, y_range::AbstractRange; ref::Union{Nothing, EquirectangularReference} = nothing)
    if ref !== nothing
        x_range, y_range = sph2xy(x_range, y_range, EQR_DEFAULT)
    end

    return norm([x_range[1], y_range[1]] - [x_range[2], y_range[2]])
end

"""
    ΔL(dist::SargassumDistribution)

Compute a spring length from a `SargassumDistribution`. This is the distance between the centers of 
diagonally-adjacent gridpoints.
"""
function ΔL(dist::SargassumDistribution)
    p1 = sph2xy(dist.lon[1], dist.lat[1], EQR_DEFAULT)
    p2 = sph2xy(dist.lon[2], dist.lat[2], EQR_DEFAULT)
    return norm(p1 - p2)
end

"""
    spring_force(xy1, xy2, parameters)

Calculate the x and y components of the force on a point particle with coordinates `xy1` 
which is attached by a spring defined by `parameters` to another point particle with coordinates `xy2`.
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
    struct InitialConditions{T}

A container for the initial conditions for a raft. 

### Fields

- `tspan`: A `Tuple{Real, Real}` such that the integration is performed for `tspan[1] ≤ t ≤ tspan[2]` where `t` is \
in days since `WATER_ITP.x.time_start`.
- `ics`: A `Vector` such that `ics[1] = n_clumps` and `ics[2:2*n_clumps+1]` represents the `[x, y]` coordinates of each clump.

### Constructors 

use `InitialConditions(;tspan, ics)`.
"""
struct InitialConditions{T<:Real}
    tspan::Tuple{T,T}
    ics::Vector{T}

    function InitialConditions(;tspan::Tuple{Real, Real}, ics::Vector{T}) where {T<:Real}
        @assert tspan[1] < tspan[2] "initial time must be less than final time"

        tspan_prom = (T(tspan[1]), T(tspan[2]))

        return new{eltype(ics)}(tspan_prom, ics)
    end
end

"""
    InitialConditions(tspan, xy0; ref = nothing)

Construct initial conditions suitable for use in `RaftParameters.ics` from a list of coordinates `xy0` of the form 
`[x1, y1, x2, y2 ..., xN, yN]`. These should be equirectangular coordinates; if `ref` is provided, the coordinates 
are converted from spherical coordinates.

Can be applied as `InitialConditions(tspan, x_range, y_range; ref = nothing)` to generate clumps in a rectangular arrangement.

Can be applied as `InitialConditions(tspan, x0, y0; ref = nothing)` for a single clump with coordinates `(x0, y0)`.
"""
function InitialConditions(tspan::Tuple{Real, Real}, xy0::Vector{<:Real}; ref::Union{Nothing, EquirectangularReference} = nothing)
    if ref !== nothing
        ics = sph2xy(xy0, ref)
    else
        ics = deepcopy(xy0)
    end

    pushfirst!(ics, length(xy0)/2)
    return InitialConditions(tspan = tspan, ics = ics)
end

function InitialConditions(
    tspan::Tuple{Real, Real},
    x_range::AbstractRange{T}, 
    y_range::AbstractRange{T}; 
    ref::Union{Nothing, EquirectangularReference} = nothing) where {T<:Real}

    @assert allunique(x_range) "`x_range` can not have repeated entries"
    @assert allunique(y_range) "`y_range` can not have repeated entries"

    if ref !== nothing
        ics_x, ics_y = sph2xy(x_range, y_range, ref)
    else
        ics_x, ics_y = x_range, y_range
    end

    ics = T[length(x_range)*length(y_range)]

    for x in ics_x, y in ics_y
        push!(ics, x, y)
    end

    return InitialConditions(tspan = tspan, ics = ics)
end

function InitialConditions(tspan::Tuple{Real, Real}, x0::Real, y0::Real; ref::Union{Nothing, EquirectangularReference} = nothing)
    if ref !== nothing
        ics = sph2xy(x0, y0, ref)
    else
        ics = [x0, y0]
    end

    pushfirst!(ics, 1)
    return InitialConditions(tspan = tspan, ics = ics)
end

"""
    InitialConditions(dist, number, sample_type, ref)

Construct [`InitialConditions`](@ref) from a `SargassumDistribution`.

### Arguments 
- `tspan`: A `Tuple{Real, Real}` such that the integration is performed for `tspan[1] ≤ t ≤ tspan[2]` where `t` is \
in days since `WATER_ITP.x.time_start`.
- `dist`: A `SargassumDistribution`.
- `weeks`: A `Vector{<:Integer}` giving the weeks of the month to consider. Each entry should be between 1 and 4 and appear only once.
- `number`: The number of clumps to initialize; interactive with `sample_type` and should be at least `1`.
- `sample_type`: A `String` identifying one of the methods of assigning clump locations based on the distribution.
    - `"levels"`: Boxes with nonzero Sargassum are divided into `number` levels of size `(minimum(dist.sargassum) - maximum(dist.sargassum))/number`.
    Each box gets a number of clumps equal to its level index. For example, if `number = 2`, then the smaller half of the boxes (by Sargassum content) 
    get 1 clump each and the larger half get 2 clumps each.
    - `"sample"`: A number `number` of samples are drawn from `dist`. Each sample is placed uniformly at random inside the corresponding box.
    - `"sorted"`: Boxes are filled with one clump placed uniformly at random inside them, starting from the box with the highest concentration. If `number` 
                is greater than the total number of boxes, repeat the loop starting again from the highest concentration box.
    - `"uniform"`: Exactly one clump is placed in the center of each box with nonzero concentration. Note that this ignores `number`.
- `ref`: An [`EquirectangularReference`](@ref). A `SargassumDistribution` has fields `lon` and `lat`, so this is necessary to 
                covert these to equirectangular coordinates.
"""
function InitialConditions(
    tspan::Tuple{Real, Real},
    dist::SargassumDistribution, 
    weeks::Vector{<:Integer},
    number::Integer,
    sample_type::String,
    ref::EquirectangularReference)

    @assert sample_type in ["sample", "sorted", "uniform", "levels"] "`sample_type` not recognized."
    @assert number > 0 "Must request at least one clump"
    @assert length(weeks) > 0 "At least one week must be selected."
    @assert all(map(x -> 1 <= x <= 4, weeks)) "Each entry of `weeks` must be between 1 and 4."
    @assert allunique(weeks) "Each week should only appear once."

    sarg = sum(dist.sargassum[:,:,week] for week in weeks)
    lons = dist.lon
    lats = dist.lat

    δ_x = abs(lons[2] - lons[1])
    δ_y = abs(lats[2] - lats[1])

    xy0 = Float64[]

    pts = Iterators.product(lons, lats) |> x -> reduce(vcat, collect(x)) # vector of (lon, lat)
    wts = sarg |> x -> reduce(vcat, x) # vector of sarg, matched with pts

    if sample_type == "sample"
        samples = Distributions.sample(pts, Weights(wts), number)

        for samp in samples
            push!(xy0, rand(Uniform(samp[1] - δ_x/2, samp[1] + δ_x/2)))
            push!(xy0, rand(Uniform(samp[2] - δ_y/2, samp[2] + δ_y/2)))
        end
    elseif sample_type == "sorted"
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
        idx = findall(x -> x > 0, wts)
        pts = pts[idx]

        for pt in pts
            push!(xy0, pt[1])
            push!(xy0, pt[2])  
        end
    elseif sample_type == "levels"
        idx = findall(x -> x > 0, wts)
        pts = pts[idx]
        wts = wts[idx]

        wtsmin, wtsmax = extrema(wts)

        for i = 1:length(pts)
            pt = pts[i]
            n_c = number*(wts[i] - wtsmin)/(wtsmax - wtsmin) |> x -> round(Integer, x, RoundUp)
            for _ = 1:n_c
                push!(xy0, rand(Uniform(pt[1] - δ_x/2, pt[1] + δ_x/2)))
                push!(xy0, rand(Uniform(pt[2] - δ_y/2, pt[2] + δ_y/2)))       
            end
        end           
    end

    xy0 = sph2xy(xy0, ref)
    n_clumps = length(xy0)/2
    pushfirst!(xy0, n_clumps)
    
    return InitialConditions(tspan = tspan, ics = xy0)
end

"""
    abstract type AbstractConnections

A supertype for all connections between clumps.

Every subtype of `AbstractConnections` should be mutable with a field `connections::Dict{U, Vector{U}} where {U<:Integer}` such 
that `connections[idx]` for an index `idx` is a vector of indices `[i1, i2, ...]` where a spring is connected 
between clumps `idx` and `[i1, i2, ...]`. This should be updated in-place as clumps grow and die, i.e. `connections` 
only shows the current connections and refers to vector indices, not absolute clump labels.

Every subtype of `AbstractConnections` should implement a `form_connections!(con::Connections, u)` method which
updates `con.connections` in place with the solution vector `u` and returns `nothing`.
"""
abstract type AbstractConnections end

"""
    mutable struct ConnectionsNone

A connection type such that no clumps are connected.

### Constructors

`ConnectionsNone()` creates an instance with no connections.
"""
mutable struct ConnectionsNone{U<:Integer} <: AbstractConnections
    connections::Dict{U, Vector{U}}

    function ConnectionsNone()
        return new{Int64}(Dict(0 => Int64[]))
    end
end

function form_connections!(con::ConnectionsNone, u::Vector{<:Real})
    n_clumps = Int64((length(u) - 1)/2)

    con.connections = Dict(1:n_clumps .=> [Int64[] for i = 1:n_clumps])    

    return nothing
end

"""
    mutable struct ConnectionsFull

A connection type such that every clump is connected to every other clump.

### Constructors

`ConnectionsFull()` creates an instance with no connections.
"""
mutable struct ConnectionsFull{U<:Integer} <: AbstractConnections
    connections::Dict{U, Vector{U}}  
    
    function ConnectionsFull()
        return new{Int64}(Dict(0 => Int64[]))
    end
end

function form_connections!(con::ConnectionsFull, u::Vector{<:Real})
    n_clumps = Int64((length(u) - 1)/2)

    con.connections = Dict(1:n_clumps .=> [[j for j = 1:n_clumps if j != i] for i = 1:n_clumps])    

    return nothing
end

"""
    mutable struct ConnectionsRadius{T}

A connection type such that every clump is connected to every clump within a given radius.

### Fields 

- `radius`: A distance in km such that each clump is connected to every clump whose distance is at most `radius` from it.

### Constructors

`ConnectionsRadius(radius)` creates an instance with no connections.
"""
mutable struct ConnectionsRadius{U<:Integer, T<:Real} <: AbstractConnections
    connections::Dict{U, Vector{U}}  
    radius::T

    function ConnectionsRadius(radius::Real)
        @assert radius > 0 "`radius` must be positive"

        return new{Int64, typeof(radius)}(Dict(0 => Int64[]), radius)
    end
end

function form_connections!(con::ConnectionsRadius, u::Vector{<:Real})
    n_clumps = Int64((length(u) - 1)/2)
    xy = reshape(@view(u[2:end]), 2, n_clumps)
    
    balltree = BallTree(xy)
    idx = inrange(balltree, xy, con.radius, true)
    idx = [filter(x -> x != i, idx[i]) for i = 1:length(idx)]
    con.connections = Dict(1:n_clumps .=> idx) 

    return nothing
end

"""
    mutable struct ConnectionsNearest{T}

A connection type such that every clump is connected to a number of its nearest neighbors.

### Fields 

- `neighbors`: The number of nearest neighbors each clump should be connected to.

### Constructors

`ConnectionsNearest(neighbors)` creates an instance with no connections.
"""
mutable struct ConnectionsNearest{U<:Integer} <: AbstractConnections
    connections::Dict{U, Vector{U}}
    neighbors::U

    function ConnectionsNearest(neighbors::Integer)
        @assert neighbors >= 0 "`neighbors` must be nonnegative"

        return new{Int64}(Dict(0 => Int64[]), neighbors)
    end
end

function form_connections!(con::ConnectionsNearest, u::Vector{<:Real})
    n_clumps = Int64((length(u) - 1)/2)
    xy = reshape(@view(u[2:end]), 2, n_clumps)

    k = con.neighbors + 1 # each point is considered one of its own nearest neighbors, so have to add 1     
    k = min(n_clumps, k) # there can not be more neighbors than there are clumps

    kdtree = KDTree(xy)
    idx, dists = knn(kdtree, xy, k, true)
    idx = [filter(x -> x != i, idx[i]) for i = 1:length(idx)]
    con.connections = Dict(1:n_clumps .=> idx)

    return nothing
end

"""
    mutable struct RaftParameters{T, U, F, C, G, L}

A container for the parameters defining a raft. Each clump and spring are identical.

### Structure 

`RaftParameters` acts as the parameter container for [`Raft!`](@ref). The solution vector `u` is a vector of length `2n_clumps + 1` 
such that `u[1]` is an "amount" parameter which controls the growth and death of clumps by biophysical effects. Then, 
`u[2*i:2*i+1]` for `i = 1:n_clumps` gives the `[x, y]` coordinates of the clump in position `i`.

### Fields
- `ics`: An [`InitialConditions`](@ref).
- `clumps`: The [`ClumpParameters`](@ref) shared by each clump in the raft.
- `springs`: The [`SpringParameters`](@ref) shared by each spring joining the clumps.
- `n_clumps_tot`: An `Integer` equal to the total number of clumps that have ever existed (i.e. it is at least the number of clumps that exist at any specific time.)
- `connections`: A subtybe of [`AbstractConnections`](@ref).
- `loc2label`: A `Dict` such that `loc2label[t]` is itself a `Dict` mapping vector indices to the absolute label of the clump in that location at
the `i`th time step. For example, `loc2label[t0][j] = j` since, at the initial time `t0`, the `j`th location contains the `j`th clump. If 
clump 1 dies at some later time `t`, then `loc2label[t][1] = 2`, `loc2label[t][2] = 3` since every clump is shifted by one to the left.
- `gd_model`: A subtype of [`AbstractGrowthDeathModel`](@ref). 
- `land`:: A subtype of [`AbstractLand`](@ref).

### Constructors 

Use `RaftParameters(; tspan, ics, clumps, springs, connections, gd_model, land)`.
The quantities `n_clumps_tot` and `loc2label` are computed automatically.
"""
mutable struct RaftParameters{T<:Real, U<:Integer, F<:Function, C<:AbstractConnections, G<:AbstractGrowthDeathModel, L<:AbstractLand}
    ics::InitialConditions{T}
    clumps::ClumpParameters{T}
    springs::SpringParameters{F, T}
    n_clumps_tot::U
    connections::C
    loc2label::Dict{T, Dict{U, U}}
    gd_model::G
    land::L

    function RaftParameters(;
        ics::InitialConditions{T},
        clumps::ClumpParameters{T},
        springs::SpringParameters{F, T},
        connections::C,
        gd_model::G,
        land::L) where {T<:Real, F<:Function, C<:AbstractConnections, G<:AbstractGrowthDeathModel, L<:AbstractLand}

        n_clumps = Int64((length(ics.ics) - 1)/2)
        loc2label = Dict(ics.tspan[1] => Dict(i => i for i = 1:n_clumps))
        form_connections!(connections, ics.ics)

        return new{T, Int64, F, C, G, L}(ics, clumps, springs, n_clumps, connections, loc2label, gd_model, land)
    end
end

function Base.show(io::IO, x::RaftParameters)
    print(io, "RaftParameters[")
    show(io, Integer(x.ics.ics[1]))
    print(io, " Clumps]")
end