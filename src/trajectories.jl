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
    RaftTrajectory(sol, rp)

Construct a [`RaftTrajectory`](@ref) from a differential equation solution `sol` and [`RaftParameters`](@ref) `rp`.

### Arguments

- `sol`: The output of `solve(Raft!, args...)`.
- `rp`: The [`RaftParameters`](@ref) used in solving the [`Raft!`](@ref) model.
- `ref`: An [`EquirectangularReference`](@ref).
"""
function RaftTrajectory(sol::AbstractMatrix, rp::RaftParameters, ref::EquirectangularReference)
    data = zeros(length(unique(sol.t)), rp.n_clumps_tot, 2)
    n_clumps_t = Int64[]
    lifetimes = Dict(i => Int64[] for i = 1:rp.n_clumps_tot)

    function update_traj!(i_u, i_t)
        u = sol.u[i_u][2:end]
        nc = Integer(length(u)/2)
        push!(n_clumps_t, nc)
        nt = length(n_clumps_t)

        for j = 1:nc
            label = rp.loc2label[sol.t[i_t]][j]
            data[nt, label, :] = u[2*j-1:2*j]
            push!(lifetimes[label], nt)
        end

        return nothing
    end

    for i = 1:length(sol.t)
        t = sol.t[i]

        if i == 1
            update_traj!(i, i)
            continue
        end

        if i == length(sol.t)
            if (sol.t[i - 1] != t)
                update_traj!(i, i)
            end

            continue
        end

        if (sol.t[i - 1] != t) 
            if (sol.t[i + 1] == t)
                # starting callbacks, so use loc2label of previous time since at current time it doesn't contain all the clumps that exist now
                update_traj!(i, i - 1)
                
            else
                update_traj!(i, i)
            end
        end
    end

    # collecting trajectories
    trajectories = Dict{Int64, Trajectory{Float64}}()
    times = unique(sol.t)
    for i = 1:rp.n_clumps_tot
        trajectories[i] = Trajectory(data[lifetimes[i],i,:], times[lifetimes[i]], ref)
    end

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