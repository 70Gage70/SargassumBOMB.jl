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
    n_clumps_max = length(keys(rp.connections))
    rp.connections[n_clumps_max + 1] = rand(keys(rp.connections), 3) # UPDATE THIS WITH CONNECTION LOGIC

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

- `land_itp`: A `StaticField2DInterpolantEQR` which gives the interpolated land locations.
"""
function die_land(land_itp::StaticField2DInterpolantEQR)

    function condition(u, t, integrator)
        return any([land_itp.u(u[2*i-1], u[2*i]) == 1.0  for i = 1:Integer(length(u)/2)])
    end

    function affect!(integrator)
        u = integrator.u
        t = integrator.t
        inds = findall([land_itp.u(u[2*i-1], u[2*i]) == 1.0  for i = 1:Integer(length(u)/2)])
        inds = [inds[i] - (i - 1) for i = 1:length(inds)]
        # if we have to delete multiple clumps in one step, deleting one clump will change the indices of the others.
        # since findall is sorted, after you delete the clump indexed by inds[1], then the clumps with indices >inds[1]
        # have their index decreased by 1, and so on

        if length(inds) == Integer(length(u)/2) # all clumps will be removed, so terminate
            terminate!(integrator)
        else
            for i in inds
                deleteat!(integrator, 2*i - 1) # e.g. index i = 2, delete the 3rd component (x coord of 2nd clump)
                deleteat!(integrator, 2*i - 1) # now the y coordinate is where the x coordinate was
                kill!(integrator.p, i, t) # remove the clump and relabel its connections
            end
        end

        println("indices $inds hit shore at time $t")
    end

    return DiscreteCallback(condition, affect!)
end


function grow_test(t_grow::Vector{<:Real})
    function condition(u, t, integrator)
        return any([abs(t - tg) < 1.0 for tg in t_grow])
    end

    function affect!(integrator)
        u = integrator.u
        t = integrator.t

        resize!(integrator, length(u) + 2)
        u[end - 1] = u[end - 1 - 2] + rand()
        u[end] = u[end - 2] + rand()

        grow!(integrator.p, t)

        println("growth")
    end

    return DiscreteCallback(condition, affect!)
end