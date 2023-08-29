include("coordinates.jl")

#################################

"""
    struct Trajectory{T}

A container for the data of a single clump's trajectory. 

### Fields

- `xy`: A `Matrix` of size `N x 2` such that `xy[i,:]` gives the `[x, y]` coordinates at the clump at time `t[i]`.
- `t`: A `Vector` of length `N` giving the time values of the trajectory.

### Constructors

Apply as `Trajectory(xy, t)` where `xy` can be a matrix or vector of vectors, `t` is a vector with the same length as `xy`.

Apply as `Trajectory(xy, t, ref)` where `ref` is an [`EquirectangularReference`](@ref) to construct the trajectory such that `xy` is converted from equirectangular to lon/lat coordinates.
"""
struct Trajectory{T<:Real}
    xy::Matrix{T}
    t::Vector{T}

    function Trajectory(xy::Union{Matrix{T}, Vector{<:Vector{T}}}, t::Vector{T}) where {T<:Real}
        if xy isa Matrix
            @assert length(t) == size(xy, 1) "`t` and `xy` must have the same length."
            @assert size(xy, 2) == 2 "`xy` must have two columns."
        
            inds = unique(i -> t[i], eachindex(t))

            return new{T}(xy[inds,:], t[inds])
        else
            @assert length(t) == length(xy) "`t` and `xy` must have the same length."

            inds = unique(i -> t[i], eachindex(t))
        
            return new{T}(stack(xy, dims = 1)[inds,:], t[inds])
        end
    end   

    function Trajectory(xy::Union{Matrix{T}, Vector{<:Vector{T}}}, t::Vector{T}, ref::EquirectangularReference) where {T<:Real}
        if xy isa Matrix
            @assert length(t) == size(xy, 1) "`t` and `xy` must have the same length."
            @assert size(xy, 2) == 2 "`xy` must have two columns."
        
            inds = unique(i -> t[i], eachindex(t))

            return new{T}(xy2sph(xy[inds,:], ref), t[inds])
        else
            @assert length(t) == length(xy) "`t` and `xy` must have the same length."

            inds = unique(i -> t[i], eachindex(t))
        
            return new{T}(xy2sph(stack(xy, dims = 1)[inds,:], ref), t[inds])
        end
    end 
end

function Base.length(tr::Trajectory)
    return length(tr.t)
end

function Base.show(io::IO, x::Trajectory)
    print(io, "Trajectory[(")
    show(io, first(x.t))
    print(io, ", ")
    show(io, last(x.t))
    print(io, "), ")
    show(io, length(x.t))
    print(io, " pts]")
end

"""
    struct RaftTrajectory{N, J}

A container for the data of a every clump's trajectory in a raft, as well as its center of mass.
    
### Fields 
- `trajectories`: A `Dict` mapping clump indices to their corresponding [`Trajectory`](@ref).
- `t`: A vector of all time possible slices across the clump trajectories.
- `n_clumps`: A vector such that `n_clumps[i]` is the number of clumps that are alive at time `t[i]`.
- `com`: A [`Trajectory`](@ref) corresponding to the center of mass of the raft. 
"""
struct RaftTrajectory{U<:Integer, T<:Real}
    trajectories::Dict{U,Trajectory{T}}
    t::Vector{T}
    n_clumps::Vector{U}
    com::Trajectory{T}
end

"""
    RaftTrajectory(sol, rp, ref)

Construct a [`RaftTrajectory`](@ref) from a differential equation solution `sol` and [`RaftParameters`](@ref) `rp`.

### Arguments

- `sol`: The output of `solve(Raft!, args...)`.
- `rp`: The [`RaftParameters`](@ref) used in solving the [`Raft!`](@ref) model.
- `ref`: An [`EquirectangularReference`](@ref).
"""
function RaftTrajectory(sol::AbstractMatrix, rp::RaftParameters, ref::EquirectangularReference)
    # unique times
    times = Float64[]

    # data[time idx, clump label, xy]
    data = zeros(length(unique(sol.t)), rp.n_clumps_tot, 2)

    # the number of clumps at each (unique) time
    n_clumps_t = Int64[]

    # map from clump label to list of time indexes it's alive
    lifetimes = Dict(i => Int64[] for i = 1:rp.n_clumps_tot)

    for i = 1:length(sol.t)
        t = sol.t[i]
        u = sol.u[i][2:end]
        nc = Integer(length(u)/2)

        if !(t in times) 
            push!(times, t)
            push!(n_clumps_t, nc)
        else
            n_clumps_t[end] = max(nc, n_clumps_t[end])
        end

        i_t = length(times)

        label_time = (i != 1) && (i != length(sol.t)) && (sol.t[i+1] == t) && (t != sol.t[i-1]) ? sol.t[i-1] : t

        for loc = 1:nc
            l2l = rp.loc2label[label_time]
            if !(loc in keys(l2l)) continue end

            label = l2l[loc]
            data[i_t, label, :] = u[2*loc-1:2*loc]
            if length(lifetimes[label]) == 0 || lifetimes[label][end] != i_t
                push!(lifetimes[label], i_t)
            end
        end
    end

    # collecting trajectories
    trajectories = Dict{Int64, Trajectory{Float64}}()

    for i = 1:rp.n_clumps_tot
        trajectories[i] = Trajectory(data[lifetimes[i],i,:], times[lifetimes[i]], ref)
    end

    # return data, n_clumps_t 

    tr_com = Trajectory(sum(data[:,i,:] for i = 1:rp.n_clumps_tot) ./ n_clumps_t, times, ref)

    return RaftTrajectory(trajectories, times, n_clumps_t, tr_com)
end

function Base.show(io::IO, x::RaftTrajectory)
    print(io, "RaftTrajectory[")
    show(io, length(keys(x.trajectories)))
    print(io, " trajectories, ")
    show(io, length(x.t))
    print(io, " times]")
end