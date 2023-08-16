using LinearAlgebra: norm

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
    spring_force(xy1, xy2, parameters)

Calculate the x and y components of the force on a point particle with coordinates `xy1` which is attached by a spring defined by `parameters` to another point particle with coordinates `xy2`.
"""
function spring_force(xy1::Vector{<:Real}, xy2::Vector{<:Real}, parameters::SpringParameters)
    d = norm(xy1 - xy2)
    return parameters.k(d)*(parameters.L/d - 1)*(xy1 - xy2)
end

"""
    mutable struct RaftParameters{C, S}

A container for the parameters defining a raft. Each clump and spring are identical.

### Fields
- `xy0`: An array representing the initial coordinates of the clumps.
- `clumps`: The [`ClumpParameters`](@ref) shared by each clump in the raft.
- `springs`: The [`SpringParameters`](@ref) shared by each spring joining the clumps.
- `connections`: A `Dict` such that `connections[idx]` is a vector of indices `idx'` where a spring is connected between clumps `idx` and `idx'`. This should be updated in-place as clumps grow and die, i.e. `connections` only shows the current connections.
- `growths`: A `Dict` such that `growths[t]` gives a vector of indices of clumps which were created at time `t`. Does not include clumps "grown" by initial conditions at time `t = 0`.
- `deaths`: A `Dict` such that `deaths[t]` gives a vector of indices of clumps which were removed at time `t`.

### Growths and Deaths

The `growths` field and `deaths` field refer to indices of the clumps which exist at that particular time, and are labelled in order along the solution vector `u`. In other words, if `deaths[t1] = [3, 7]`, this means that `u[5:6]` and `u[13:14]` were both deleted at time `t1`. If `t2 > t1` and `deaths[t2] = [3, 7]`, it again means that `u[5:6]` and `u[13:14]` were both deleted at time `t2`, however, these would now correspond to different actual clumps since the deletion at time `t1` shifted the indices.
"""
mutable struct RaftParameters{
    T<:AbstractArray, 
    C<:ClumpParameters, 
    S<:SpringParameters, 
    N1<:AbstractDict,
    N2<:AbstractDict
    }

    xy0::T
    clumps::C
    springs::S
    connections::N1
    growths::N2
    deaths::N2
end

"""
    RaftParameters(x_range, y_range, clump_parameters, spring_parameters; network_type, name)

Construct [`RaftParameters`](@ref) in a rectangular arrangement.

The initial conditions `xy0` are collected in a vector of length `n_row x n_col x 2`. And are arranged as `[x1, y1, x2, y2 ...]`.

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
function RaftParameters(
    x_range::AbstractRange{<:Real}, 
    y_range::AbstractRange{<:Real},
    clump_parameters::ClumpParameters, 
    spring_k::Function; 
    network_type::String = "nearest")

    @assert network_type in ["nearest", "full", "none"] "`network_type` not recognized."
    @assert step(x_range) > 0 "x range should be inreasing."
    @assert step(y_range) > 0 "y range should be increasing."

    # a rectangular mesh
    network = reverse.(collect(Iterators.product(reverse(y_range), x_range))) # reverse so that the first row has the largest y
    network = reshape(stack(network, dims = 1), (size(network, 1), size(network, 2), 2)) # convert from matrix of tuples to array
    n_col = length(x_range)
    n_row = length(y_range)
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

    xy0 = [network[i, j, k] for i = 1:n_row for j = 1:n_col for k = 1:2]

    n(i, j) = (i - 1) * n_col + j

    connections_flat = Dict{Int64, Vector{Int64}}()

    for key in keys(connections)
        connections_flat[n(key...)] = connections[key] .|> x -> n(x...)
    end

    growths = Dict{Float64, Vector{Int64}}()
    deaths = Dict{Float64, Vector{Int64}}()    

    return RaftParameters(xy0, clump_parameters, spring_parameters, connections_flat, growths, deaths)
end


"""
    kill!(rp::RaftParameters, i, t)

Remove the clump with index `i` and its connections from `rp` and update `rp.deaths` at time `t.`

Indices whose value is greater than i are then shifted down by 1.
"""
function kill!(rp::RaftParameters, i::Integer, t::Float64)
    delete!(rp.connections, i) # remove i from keys
    rp.connections = Dict(a => filter(x -> x != i, b) for (a,b) in rp.connections) # remove i from values

    less_i(x) = x > i ? x - 1 : x
    rp.connections = Dict(less_i(a) => less_i.(b) for (a,b) in rp.connections) # shift every label >i down by 1

    if t in keys(rp.deaths)
        push!(rp.deaths[t], i)
    else
        rp.deaths[t] = [i]
    end

    return nothing
end

"""
    grow!(rp::RaftParameters, t)

Blah blah blah.
"""
function grow!(rp::RaftParameters, t::Float64)
    n_clumps_max = maximum(keys(rp.connections))
    rp.connections[n_clumps_max + 1] = [] # UPDATE THIS WITH CONNECTION LOGIC

    if t in keys(rp.growths)
        push!(rp.growths[t], n_clumps_max + 1)
    else
        rp.growths[t] = [n_clumps_max + 1]
    end
    
    return nothing
end

"""
    raft_trajectories(sol, rp)

Construct a `Dict`, `rt`, such that `rt[i]` is a `N x 3` matrix giving the trajectory of the `i`th clump.

The convention is `rt[i][x, y, t]`.
"""
function raft_trajectories(sol::AbstractMatrix, rp::RaftParameters)
    n_total_clumps = Integer(length(sol[1])/2) # keeps track of the total number of clumps that have ever existed

    # loc_to_label keeps track of which clump actually has its coordinates in positions u[2*i-1, 2*i].
    # when initialized, this is exactly clump i, but it will change after growths and deaths
    loc_to_label = [i for i = 1:n_total_clumps]

    # placing and labeling initial clumps
    tr = Dict{Int64, Matrix{Float64}}()
    for i = 1:Integer(length(sol[1])/2)
        tr[i] = [sol[1][2*i-1] sol[1][2*i] sol.t[1]]
    end

    for j = 2:length(sol)
        if sol.t[j] != sol.t[j - 1] # there was no growth or death at this time
            for i = 1:Integer(length(sol[j])/2)
                if loc_to_label[i] in keys(tr) # clump trajectory is already started, so add to it
                    tr[loc_to_label[i]] = vcat(tr[loc_to_label[i]], [sol[j][2*i-1] sol[j][2*i] sol.t[j]])
                else # clump trajectory must be started
                    tr[loc_to_label[i]] = [
                                            sol[j-1][2*i-1] sol[j-1][2*i] sol.t[j-1]; # need previous time since clump was "born" then
                                            sol[j][2*i-1] sol[j][2*i] sol.t[j]]
                end
            end 
        else # note that a growth AND a death could happen in the same step; handle deaths first
            if sol.t[j] in keys(rp.deaths)
                # suppose that rp.deaths[sol.t[j]] = [k], 
                # then loc_to_label[k'] should be set to loc_to_label[k' +  1] for each  k <= k <= end-1 (i.e. move everything to the left)
                # and the last entry should be removed
                for k1 in rp.deaths[sol.t[j]]
                    for k2 = k1:length(loc_to_label)-1
                        loc_to_label[k2] = loc_to_label[k2 + 1]
                    end
                    
                    pop!(loc_to_label)
                end
            end

            if sol.t[j] in keys(rp.growths)
                # add an extra entry at the end of loc_to_label and increment n_total_clumps for each growth
                for k1 in rp.growths[sol.t[j]]
                    n_total_clumps = n_total_clumps + 1
                    push!(loc_to_label, n_total_clumps)
                end
            end
        end
    end

    return tr
end