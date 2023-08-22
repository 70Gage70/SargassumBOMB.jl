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
    # n_total_clumps keeps track of the total number of clumps that have ever existed
    # increases with growths and doesn't change with deaths
    n_total_clumps = Integer(length(sol[1])/2) 

    # n_clumps keeps track of the number of clumps at each time slice
    # n_clumps[1] = the initial number of clumps = n_total_clumps
    n_clumps = [n_total_clumps]

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

            push!(n_clumps, Integer(length(sol[j])/2))

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

    # collecting everything into a trajectory
    trajectories = Dict{Int64, Trajectory{Float64}}()
    for i in keys(tr)
        trajectories[i] = Trajectory(tr[i][:,1:2], tr[i][:,3], ref)
    end
    
    # all unique times
    times = unique(sol.t)

    # arrange trajectories into one big array such that they are matched up by time slice
    # this makes computing the COM much easier
    com_array = zeros(n_total_clumps, length(times), 2)

    for i = 1:n_total_clumps
        # first(trajectories[i].t) is the first time this clump exists
        # `ind` is therefore the index where this trajectory should be slotted into `com_array`
        ind = searchsortedfirst(times, first(trajectories[i].t)) 
        com_array[i,ind:ind+length(trajectories[i])-1,:] .= trajectories[i].xy
    end

    com = [sum(com_array[ci,ti,:] for ci = 1:n_total_clumps)/n_clumps[ti] for ti = 1:length(times)]

    return RaftTrajectory(trajectories, times, n_clumps, Trajectory(com, times))
end

function Base.show(io::IO, x::RaftTrajectory)
    print(io, "RaftTrajectory[")
    show(io, length(keys(x.trajectories)))
    print(io, " trajectories, ")
    show(io, length(x.t))
    print(io, " times]")
end