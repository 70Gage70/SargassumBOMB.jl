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
    cb_land(land::AbstractLand)

Create a `DiscreteCallback` which updates `integrator.p` at the end of each time step using the 
[`AbstractLand`](@ref) in `land`.
"""
function cb_land(land::AbstractLand)
    return DiscreteCallback(land, land)
end

"""
    cb_growth_death(model::AbstractGrowthDeathModel)

Create a `DiscreteCallback` which updates `integrator.p` at the end of each time step using the 
[`AbstractGrowthDeathModel`](@ref) in `model`.
"""
function cb_growth_death(model::AbstractGrowthDeathModel)
    return DiscreteCallback(model, model)
end

"""
    cb_connections()

Create a `DiscreteCallback` which updates `integrator.p.connections` at the end of each time step using its `form_connections!` method.

This should be the last callback in a `CallbackSet`.
"""
function cb_connections()
    condition!(u, t, integrator) = true    

    function affect!(integrator)
        form_connections!(integrator.p.connections, integrator.u)
        return nothing
    end

    return DiscreteCallback(condition!, affect!, save_positions = (false, false))
end

"""
    kill!(integrator, i)

Remove the clump with index `i` (by vector location) from `integrator.u`. 

For the [`RaftParameters`](@ref), `rp = integrator.p` update `rp.loc2label` and `rp.gd_model` appropriately.
"""
function kill!(integrator::SciMLBase.DEIntegrator, i::Integer)
    # if this is the last clump, terminate the integration
    if n_clumps(integrator.u) == 1
        terminate!(integrator)
        return nothing
    end

    # first remove the appropriate u elements from `integrator`
    deleteat!(integrator, 2*i-1) # e.g. index i = 2, delete the 3rd component (x coord of 2nd clump)
    deleteat!(integrator, 2*i-1) # now the y coordinate is where the x coordinate was

    # update rp.loc2label at the current time
    rp = integrator.p 
    t = integrator.t
    gtr_i(x) = x >= i ? x + 1 : x
    rp.loc2label[t] = Dict(j => rp.loc2label[t][gtr_i(j)] for j = 1:n_clumps(integrator.u))

    # update rp.gd_model.S at the current time
    deleteat!(rp.gd_model.S, i)

    return nothing
end

"""
    kill!(integrator, inds)

Call [`kill!(integrator, i)`](@ref) on each element of `inds`.    

This removes several clumps "simultaneously" while taking into account the fact that removing a clump shifts the clump to which each element of `inds` refers.
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
    grow!(integrator, location)

Add a clump to the [`RaftParameters`](@ref), `rp = integrator.p` with an index equal to `rp.n_clumps_tot + 1` and also update `rp.loc2label` appropriately.

### Location 

`location` can be a pre-defined flag, an integer, or a `[x, y]` vector. The default value is the flag `"parent"`.

The possible flags are:
- `"parent"`: A parent clump is chosen randomly among clumps that already exist, and the new clump is placed a distance `integrator.rp.springs.L` away and at a 
random angle from it.
- `"com"`: The same as `"parent"`, except the centre location is at the center of mass of the raft.

If `location` is an `Integer` with value `i`, then the new clump will be grown with `i`th clump (by vector location) as its parent.

If `location` is a `Vector{<:Real}`, the new clump will be placed at those `[x, y]` coordinates. 
"""
function grow!(
    integrator::SciMLBase.DEIntegrator; 
    location::Union{String, Integer, Vector{<:Real}} = "parent")
    
    rp = integrator.p
    u = integrator.u
    n_clumps_old = n_clumps(u)
    n_clumps_new = n_clumps_old + 1

    # resize the integrator and add the new (x, y) components
    resize!(integrator, length(u) + 2)

    if location isa String 
        @argcheck location in ["parent", "com"] "If `location` is a string, it must be either $("parent") or $("com")."
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
        @argcheck length(location) == 2 "The location vector must be [x, y] coordinates."
        u[end-1:end] .= location
    end

    # update loc2label and n_clumps_tot
    t = integrator.t
    rp.n_clumps_tot = rp.n_clumps_tot + 1
    rp.loc2label[t][n_clumps_new] = rp.n_clumps_tot 

    # update rp.gd_model.S
    push!(rp.gd_model.S, 0.0)
    
    return nothing
end
