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
    n_c = n_clumps(u)

    con.connections = Dict(1:n_c .=> [Int64[] for i = 1:n_c])    

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
    n_c = n_clumps(u)

    con.connections = Dict(1:n_c .=> [[j for j = 1:n_c if j != i] for i = 1:n_c])    

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
        @argcheck radius > 0 "`radius` must be positive"

        return new{Int64, typeof(radius)}(Dict(0 => Int64[]), radius)
    end
end

function form_connections!(con::ConnectionsRadius, u::Vector{<:Real})
    n_c = n_clumps(u)
    xy = reshape(view(u,:), 2, n_c)
    
    balltree = BallTree(xy)
    idx = inrange(balltree, xy, con.radius, true)
    idx = [filter(x -> x != i, idx[i]) for i = 1:length(idx)]
    con.connections = Dict(1:n_c .=> idx) 

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
        @argcheck neighbors >= 0 "`neighbors` must be nonnegative"

        return new{Int64}(Dict(0 => Int64[]), neighbors)
    end
end

function form_connections!(con::ConnectionsNearest, u::Vector{<:Real})
    n_c = n_clumps(u)
    xy = reshape(view(u,:), 2, n_c)

    k = con.neighbors + 1 # each point is considered one of its own nearest neighbors, so have to add 1     
    k = min(n_c, k) # there can not be more neighbors than there are clumps

    kdtree = KDTree(xy)
    idx, dists = knn(kdtree, xy, k, true)
    idx = [filter(x -> x != i, idx[i]) for i = 1:length(idx)]
    con.connections = Dict(1:n_c .=> idx)

    return nothing
end


"""
    abstract type AbstractSpring

A supertype for all spring parameters. Each clump, when conncted, is joined by the same kind of spring.

Every subtype of `AbstractSpring` should have a field `k::Function` representing the stiffness force 
and callable as `k(x)` as well as a field `L::Real` representing the spring's natural length.

All forces are computed using `parameters.k(d)*(parameters.L/d - 1)*(xy1 - xy2)` where `d = norm(xy1 - xy2)`.
"""
abstract type AbstractSpring end

"""
    HookeSpring{T, Tk}

A subtype of `AbstractSpring` representing a spring with a constant stiffness.

### Constructor

`HookeSpring(k::Real, L::Real)`
"""
struct HookeSpring{T<:Real, Tk<:Function} <: AbstractSpring
    k::Tk
    L::T

    function HookeSpring(k::Real, L::Real)
        sk(x::Real; k::Real = k) = k
        return new{typeof(L), typeof(sk)}(sk, L)
    end
end

"""
    BOMBSpring{T, Tk}

A subtype of `AbstractSpring` representing a BOMB spring of the form `A * (exp((x - 2*L)/0.2) + 1)^(-1)`.

### Extra fields

- `A`: The amplitude of the force.

### Constructor

`BOMBSpring(A::Real, L::Real)`
"""
struct BOMBSpring{T<:Real, Tk<:Function} <: AbstractSpring
    k::Tk
    L::T
    A::T

    function BOMBSpring(A::Real, L::Real)
        sk(x::Real; A::Real = A, L::Real = L) = A * (exp((x - 2*L)/0.2) + 1)^(-1)
        return new{typeof(L), typeof(sk)}(sk, L, A)
    end
end

"""
    ΔL(x_range, y_range; to_xy)

Compute a spring length from a rectangular arrangement of clumps provided by `x_range` and `y_range`. This is the distance between the centers of 
diagonally-adjacent gridpoints.

These should be equirectangular coordinates; if `to_xy == true` the ranges are 
converted from spherical to equirectangular coordinates. Default `false`.
"""
function ΔL(x_range::AbstractRange, y_range::AbstractRange; to_xy::Bool = false)
    if to_xy
        x_range, y_range = sph2xy(x_range, y_range)
    end

    return norm([x_range[1], y_range[1]] - [x_range[2], y_range[2]])
end

"""
    ΔL(dist::SargassumDistribution)

Compute a spring length from a `SargassumDistribution`. This is the equirectangular distance between the centers of 
diagonally-adjacent gridpoints.
"""
function ΔL(dist::SargassumDistribution)
    p1 = sph2xy(dist.lon[1], dist.lat[1])
    p2 = sph2xy(dist.lon[2], dist.lat[2])
    return norm(p1 - p2)
end

"""
    ΔL(ics::InitialConditions)

Compute a spring length from a `InitialConditions`. This is the median among all pairwise 
equirectangular distances between points' 5 nearest neighbors.
"""
function ΔL(ics::InitialConditions)
    xy = reshape(view(ics.ics,:), 2, n_clumps(ics.ics))

    kdtree = KDTree(xy)
    k = 5
    idx, dists = knn(kdtree, xy, k, true)
    dists = [sum(d)/(k - 1) for d in dists] # k - 1 since self is included   

    return median(dists)
end