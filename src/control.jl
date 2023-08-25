using StatsBase
using Distributions

############################################################

"""
    n_clumps(u)

Return the number of clumps in the solution vector `u`. This is `Integer(length(u)/2)`.
"""
function n_clumps(u::Vector{<:Real})
    return Integer(length(u)/2)
end

"""
    n_clumps(u)

Return the `[x, y]` coordinates of the `i`th clump in the solution vector `u`. This is `u[2*i - 1:2*i]`.
"""
function clump_i(u::Vector{<:Real}, i::Integer)
    return u[2*i - 1:2*i]
end

"""
    com(u)

Return the center of mass `[x, y]` coordinates of the solution vector `u`.
"""
function com(u::Vector{<:Real})
    return [mean(u[1:2:end]), mean(u[2:2:end])]
end

"""
    kill!(integrator, i)

Remove the clump with index `i` from `integrator.u`. Remap the connections from the [`RaftParameters`](@ref), `rp = integrator.p` and update `rp.deaths` at time `integrator.t.`

Indices whose value is greater than i are then shifted down by 1.
"""
function kill!(integrator::SciMLBase.DEIntegrator, i::Integer)
    # if this is the last clump, terminate the integration
    if length(integrator.u) == 2
        terminate!(integrator)
        return nothing
    end

    # first remove the appropriate u elements from `integrator`
    deleteat!(integrator, 2*i - 1) # e.g. index i = 2, delete the 3rd component (x coord of 2nd clump)
    deleteat!(integrator, 2*i - 1) # now the y coordinate is where the x coordinate was

    # now remap the connections in RaftParameters
    rp = integrator.p
    
    delete!(rp.connections, i) # remove i from keys
    rp.connections = Dict(a => filter(x -> x != i, b) for (a,b) in rp.connections) # remove i from values

    less_i(x) = x > i ? x - 1 : x
    rp.connections = Dict(less_i(a) => less_i.(b) for (a,b) in rp.connections) # shift every label >i down by 1

    # finally, update rp.deaths at the current time
    t = integrator.t

    if t in keys(rp.deaths)
        push!(rp.deaths[t], i)
    else
        rp.deaths[t] = [i]
    end

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
    grow!(integrator, x, y, connections)

Add a clump to the [`RaftParameters`](@ref), `rp = integrator.p` with an index equal to the maximum clump index and update `rp.growths` at time `intergrator.t`.

### Locations 

`locations` can be a pre-defined flag, and integer, or a `[x, y]` vector.

The possible flags are:
- `"parent"`: A parent clump is chosen randomly among clumps that already exist, and the new clump is placed a distance `integrator.rp.springs.L` away and at a 
random angle from it.
- `"com"`: The same as `"parent"`, except the centre location is at the center of mass of the raft.

If `locations` is an `Integer` with value `i`, then the new clump will be grown with `i`th clump as its parent.

If `locations` is a `Vector{<:Real}`, the new clump will be placed at those `[x, y]` coordinates. 

### Connections 

`connections` can be a pre-defined flag, an integer, or a vector of integers.

The possible flags are:
- `"full"`: The new clump is connected to each other clump.
- `"none"`: The new clump it not connected to any other clump.

If `connections` is an `Integer` with value `n`, the new clump will be connected to the nearest `n` clumps. If `n` is ≥ the total number of clumps, 
this is equivalent to `"full"`.

If `connections` is a `Vector{<:Integer}`, the new clump will be connected to clumps with those indices.
"""
function grow!(
    integrator::SciMLBase.DEIntegrator, 
    locations::Union{String, Integer, Vector{<:Real}}, 
    connections::Union{String, Integer, Vector{<:Integer}})
    
    rp = integrator.p
    u = integrator.u

    # resize the integrator and add the new components
    resize!(integrator, length(u) + 2)

    if locations isa String 
        @assert locations in ["parent", "com"] "If `connections` is a string, it must be either $("parent") or $("com")."
        if locations == "parent"
            parent = rand(keys(rp.connections))
            r, θ = rp.springs.L, rand(Uniform(0, 2*π))
            u[end-1:end] = clump_i(u, parent) + [r*cos(θ), r*sin(θ)]
        elseif locations == "com"
            r, θ = rp.springs.L, rand(Uniform(0, 2*π))
            u[end-1:end] = com(u) + [r*cos(θ), r*sin(θ)]
        end
    elseif locations isa Integer
        parent = locations
        r, θ = rp.springs.L, rand(Uniform(0, 2*π))
        u[end-1:end] = clump_i(u, parent) + [r*cos(θ), r*sin(θ)]
    elseif locations isa Vector 
        @assert length(locations) == 2 "The locations vector must be [x, y] coordinates."
        u[end-1:end] .= locations
    end

    # update the connections, this new clump's label should be the highest current label + 1
    n_clumps_max = length(keys(rp.connections))

    if connections isa String 
        @assert connections in ["full", "none"] "If `connections` is a string, it must be either $("full") or $("none")."
        if connections == "full"
            rp.connections[n_clumps_max + 1] = collect(1:n_clumps_max)
        elseif connections == "none"
            rp.connections[n_clumps_max + 1] = valtype(rp.connections)()
        end
    elseif connections isa Integer
        dists = [norm([x, y] - clump_i(u, i)) for i = 1:n_clumps_max-1]
        closest = partialsortperm(dists, 1:min(connections, n_clumps_max), rev = true) |> collect
        rp.connections[n_clumps_max + 1] = closest
    elseif connections isa Vector
        rp.connections[n_clumps_max + 1] = closest
    end

    # update rp.growths at the current time
    t = integrator.t
    if t in keys(rp.growths)
        push!(rp.growths[t], n_clumps_max + 1)
    else
        rp.growths[t] = [n_clumps_max + 1]
    end
    
    return nothing
end


"""
    growth_death_temperature

IN DEVELOPMENT (SLOW).

Checks if the average temperature in the last day was above or below T_opt and whether there has been a T-growth/death in the past day.
"""
function growth_death_temperature(temp_itp::InterpolatedField; t_lag::Real, T_opt::Real, T_thresh::Real)
    function condition(u, t, integrator)
        rp = integrator.p

        if t <= t_lag
            return false
        end

        xy_now = integrator.sol(t)
        xy_past = integrator.sol(t - t_lag)

        com_now = [mean(xy_now[1:2:end]), mean(xy_now[2:2:end])]
        com_past = [mean(xy_past[1:2:end]), mean(xy_past[2:2:end])]

        temp = (temp_itp.u(com_now..., t) + temp_itp.u(com_past..., t-t_lag))/2

        recent_growth = t - maximum(keys(rp.growths), init = 0.0) > t_lag
        recent_death = t - maximum(keys(rp.deaths), init = 0.0) > t_lag
        critical_temp = abs(temp - T_opt) > T_thresh

        return recent_growth && recent_death && critical_temp
    end
    
    function affect!(integrator)
        rp = integrator.p
        t = integrator.t

        xy_now = integrator.sol(t)
        xy_past = integrator.sol(t - t_lag)

        com_now = [mean(xy_now[1:2:end]), mean(xy_now[2:2:end])]
        com_past = [mean(xy_past[1:2:end]), mean(xy_past[2:2:end])]

        temp = (temp_itp.fields[:temp](com_now..., t) + temp_itp.fields[:temp](com_past..., t-t_lag))/2

        if temp < T_opt
            kill_ind = rand(keys(rp.connections))
            kill!(integrator, kill_ind)

            @info "Clump $kill_ind T-death at $t"

        elseif temp > T_opt
            grow!(integrator)

            @info "T-growth at $t"
        end

        return nothing
    end

    return DiscreteCallback(condition, affect!)
end

function grow_test(t_grow::Vector{<:Real})
    function condition(u, t, integrator)
        return any([abs(t - tg) < 1.0 for tg in t_grow])
    end

    function affect!(integrator)
        grow!(integrator)

        println("growth")
    end

    return DiscreteCallback(condition, affect!)
end




