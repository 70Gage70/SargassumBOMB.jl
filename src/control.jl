using StatsBase
using Distributions

############################################################

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
    grow!(integrator)

Add a clump to the [`RaftParameters`](@ref), `rp = integrator.p` with an index equal to the maximum clump index and update `rp.growths` at time `intergrator.t`.

### Growth logic

To grow a new clump, a currently existing clump is chosen uniformly at random, then the new clump is placed a distance `rp.springs.L` away from that clump in a random direction and attached with one spring.
"""
function grow!(integrator::SciMLBase.DEIntegrator)
    rp = integrator.p
    u = integrator.u
    
    # first determine the location of the new clump
    mother = rand(keys(rp.connections))
    x, y = u[2*mother-1:2*mother]
    r, θ = rp.springs.L, rand(Uniform(0, 2*π))
    x, y = x + r*cos(θ), y + r*sin(θ)

    # resize the integrator and add the new components
    resize!(integrator, length(u) + 2)
    u[end-1:end] .= x, y

    # update the connections, this new clump's label should be the highest current label + 1
    n_clumps_max = length(keys(rp.connections))
    rp.connections[n_clumps_max + 1] = [mother]

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
    die_land(land_itp)

Create a `DiscreteCallback` which kills clumps when they reach the shore and terminates the integration if no clumps remain.

### Arguments

- `land_itp`: A `InterpolatedField` which gives the interpolated land locations.
"""
function die_land(land_itp::InterpolatedField)

    function condition(u, t, integrator)
        return any([land_itp.fields[:land](u[2*i-1], u[2*i]) == 1.0  for i = 1:Integer(length(u)/2)])
    end

    function affect!(integrator)
        xy = integrator.u
        inds = findall([land_itp.fields[:land](xy[2*i-1], xy[2*i]) == 1.0  for i = 1:Integer(length(xy)/2)])
        kill!(integrator, inds)

        @info "Clump $inds hit shore at time $(integrator.t)"
    end

    return DiscreteCallback(condition, affect!)
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


function n_clumps(u::Vector{<:Vector{<:Real}})
    return Integer(length(u)/2)
end

function n_clumps(integrator::SciMLBase.DEIntegrator)
    return Integer(length(integrator.u)/2)
end

function clump_i(u::Vector{<:Vector{<:Real}}, i::Integer)
    return u[2*i - 1:2*i]
end

function clump_i(integrator::SciMLBase.DEIntegrator, i::Integer)
    return integrator.u[2*i - 1:2*i]
end

