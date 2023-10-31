using StatsBase
using Distributions
using NearestNeighbors
using Crayons.Box

############################################################

"""
    n_clumps(u)

Return the number of clumps in the solution vector `u`. This is `floor(Int64, length(u)/2)`.
"""
function n_clumps(u::Vector{<:Real})
    return floor(Int64, length(u)/2)
end

"""
    clump_i(u, i)

Return the `[x, y]` coordinates of the `i`th clump in the solution vector `u`. This is `u[2*i:2*i+1]`.
"""
function clump_i(u::Vector{<:Real}, i::Integer)
    return u[2*i:2*i+1]
end

"""
    com(u)

Return the center of mass `[x, y]` coordinates of the solution vector `u`.
"""
function com(u::Vector{<:Real})
    return [mean(u[2:2:end]), mean(u[3:2:end])]
end

"""
    cb_update(; showprogress = false)

This callback is mandatory, and must be the first callback in a `CallbackSet`.

Create a `DiscreteCallback` which updates `integrator.p.loc2label` at the end of each time step using a `deepcopy` of the previous step.

If `showprogress = true`, then the percentage completion of the integration will be displayed.


"""
function cb_update(;showprogress::Bool = false)
    condition!(u, t, integrator) = true   

    function affect!(integrator)
        if showprogress
            t0, tend = integrator.sol.prob.tspan
            val = round(100*(integrator.t - t0)/(tend - t0), sigdigits = 3)
            print(WHITE_BG("Integrating: $(val)%   \r"))
            flush(stdout)
        end

        integrator.p.loc2label[integrator.t] = deepcopy(integrator.p.loc2label[integrator.tprev])
        return nothing
    end

    return DiscreteCallback(condition!, affect!, save_positions = (false, false))
end

"""
    cb_connections(;network_type, neighbor_parameter)

Create a `DiscreteCallback` which updates `integrator.p.connections` at the end of each time step using [`form_connections`](@ref).

The optional arguments are identical to those of [`form_connections`](@ref). THe default behavior is for each clump to 
establish connections with the 10 closest clumps.

This should be the last callback in a `CallbackSet`.
"""
function cb_connections(;network_type::String = "nearest", neighbor_parameter::Union{Nothing, Real} = 10)
    condition!(u, t, integrator) = true    

    function affect!(integrator)
        integrator.p.connections = form_connections(integrator.u, network_type, neighbor_parameter = neighbor_parameter)
        return nothing
    end

    return DiscreteCallback(condition!, affect!, save_positions = (false, false))
end

"""
    kill!(integrator, i)

Remove the clump with index `i` from `integrator.u`. For the [`RaftParameters`](@ref), `rp = integrator.p` update `rp.connections` and `rp.loc2label` appropriately.
"""
function kill!(integrator::SciMLBase.DEIntegrator, i::Integer)
    # if this is the last clump, terminate the integration
    if length(integrator.u) == 3
        terminate!(integrator)
        return nothing
    end

    # first remove the appropriate u elements from `integrator`
    deleteat!(integrator, 2*i) # e.g. index i = 2, delete the 4th component (x coord of 2nd clump)
    deleteat!(integrator, 2*i) # now the y coordinate is where the x coordinate was

    # now remap the connections in RaftParameters
    rp = integrator.p
    
    delete!(rp.connections, i) # remove i from keys
    rp.connections = Dict(a => filter(x -> x != i, b) for (a,b) in rp.connections) # remove i from values

    less_i(x) = x > i ? x - 1 : x
    rp.connections = Dict(less_i(a) => less_i.(b) for (a,b) in rp.connections) # shift every label >i down by 1

    # update rp.loc2label at the current time 
    t = integrator.t
    gtr_i(x) = x >= i ? x + 1 : x
    rp.loc2label[t] = Dict(j => rp.loc2label[t][gtr_i(j)] for j = 1:n_clumps(integrator.u))

    # update u[1] to remove a clump
    integrator.u[1] = integrator.u[1] - 1

    return nothing
end

"""
    kill!(integrator, inds)

Call [`kill!(integrator, i)`](@ref) on each element of `inds`.    

This removes several clumps "simultaneously" while taking into account the fact that removing a clump shifts the clump to which each element of `inds` references.
"""
function kill!(integrator::SciMLBase.DEIntegrator, inds::Vector{<:Integer})
    if length(inds) == 0
        return nothing
    end

    # if we have to delete multiple clumps in one step, deleting one clump will change the indices of the others.
    # that is, after you delete the clump indexed by inds[1], then the clumps with indices >inds[1]
    # have their index decreased by 1, and so on
    sort!(inds)    
    inds = [inds[i] - (i - 1) for i = 1:length(inds)]

    for i in inds
        kill!(integrator, i)
    end

    return nothing
end


"""
    grow!(integrator, location, connections)

Add a clump to the [`RaftParameters`](@ref), `rp = integrator.p` with an index equal to `rp.n_clumps_tot + 1` and also update `rp.connections` and `rp.loc2label` appropriately.

### Location 

`location` can be a pre-defined flag, an integer, or a `[x, y]` vector.

The possible flags are:
- `"parent"`: A parent clump is chosen randomly among clumps that already exist, and the new clump is placed a distance `integrator.rp.springs.L` away and at a 
random angle from it.
- `"com"`: The same as `"parent"`, except the centre location is at the center of mass of the raft.

If `location` is an `Integer` with value `i`, then the new clump will be grown with `i`th clump (by vector location) as its parent.

If `location` is a `Vector{<:Real}`, the new clump will be placed at those `[x, y]` coordinates. 

### Connections 

`connections` can be a pre-defined flag, an integer, or a vector of integers.

The possible flags are:
- `"full"`: The new clump is connected to each other clump.
- `"none"`: The new clump is not connected to any other clump.

If `connections` is an `Integer` with value `n`, the new clump will be connected to the nearest `n` clumps. If `n` is ≥ the total number of clumps, 
this is equivalent to `"full"`.

If `connections` is a `Vector{<:Integer}`, the new clump will be connected to clumps with those indices (by vector location).
"""
function grow!(
    integrator::SciMLBase.DEIntegrator; 
    location::Union{String, Integer, Vector{<:Real}}, 
    connections::Union{String, Integer, Vector{<:Integer}})
    
    rp = integrator.p
    u = integrator.u
    n_clumps_old = n_clumps(u)
    n_clumps_new = n_clumps_old + 1

    # resize the integrator and add the new components
    resize!(integrator, length(u) + 2)

    if location isa String 
        @assert location in ["parent", "com"] "If `connections` is a string, it must be either $("parent") or $("com")."
        if location == "parent"
            parent = rand(1:n_clumps_old)
            r, θ = rp.springs.L, rand(Uniform(0, 2*π))
            u[end-1:end] = clump_i(u, parent) + [r*cos(θ), r*sin(θ)]
        elseif location == "com"
            r, θ = rp.springs.L, rand(Uniform(0, 2*π))
            u[end-1:end] = com(u) + [r*cos(θ), r*sin(θ)]
        end
    elseif location isa Integer
        parent = location
        r, θ = rp.springs.L, rand(Uniform(0, 2*π))
        u[end-1:end] = clump_i(u, parent) + [r*cos(θ), r*sin(θ)]
    elseif location isa Vector 
        @assert length(location) == 2 "The location vector must be [x, y] coordinates."
        u[end-1:end] .= location
    end

    # update the connections
    if connections isa String 
        @assert connections in ["full", "none"] "If `connections` is a string, it must be either $("full") or $("none")."
        if connections == "full"
            rp.connections[n_clumps_new] = collect(1:n_clumps_old)
        elseif connections == "none"
            rp.connections[n_clumps_new] = valtype(rp.connections)()
        end
    elseif connections isa Integer
        dists = [norm([x, y] - clump_i(u, i)) for i = 1:n_clumps_old-1]
        closest = partialsortperm(dists, 1:min(connections, n_clumps_old), rev = true) |> collect
        rp.connections[n_clumps_new] = closest
    elseif connections isa Vector
        rp.connections[n_clumps_new] = connections
    end

    # update loc2label and n_clumps_tot
    t = integrator.t
    rp.n_clumps_tot = rp.n_clumps_tot + 1
    rp.loc2label[t][n_clumps_new] = rp.n_clumps_tot 

    # update u[1] with the new clump
    integrator.u[1] = integrator.u[1] + 1
    
    return nothing
end

function grow_test(t_grow::Vector{<:Real})
    function condition(u, t, integrator)
        return any([abs(t - tg) < 1.0 for tg in t_grow])
    end

    function affect!(integrator)
        grow!(integrator, location = "parent", connections = "full")

        println("Growing $([integrator.p.n_clumps_tot]) at time $(integrator.t)")
    end

    return DiscreteCallback(condition, affect!)
end




