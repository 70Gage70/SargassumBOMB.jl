"""
    abstract type AbstractGrowthDeathModel

The abstract type for growth and death models.
    
Subtypes must have a field `S`, a vector of length `n_clumps(u)` representing an "amount" or "mass" for each clump.

Subtypes must have a field `S0`, a vector of possible initialization values for `S`. For example, if `S0 = [1.0]`, then
each new clump will start with `S[i] = 1.0`. If `S0 = [1.0, 2.0]` then each new clump will start with one of `1.0` or 
`2.0` chosen uniformly at random.
"""
abstract type AbstractGrowthDeathModel end 

"""
    struct ImmortalModel{T}

An `AbstractGrowthDeathModel` such that no growth or death occurs.

### Constructors 

Use `ImmortalModel(ics::InitialConditions; S0 = [0.0]).`
"""
struct ImmortalModel{T<:Real} <: AbstractGrowthDeathModel
    S::Vector{T}
    S0::Vector{T}

    function ImmortalModel(ics::I; S0::Vector{T} = [0.0]) where {I<:InitialConditions, T<:Real}
        S = fill(rand(S0), length(ics.ics))
        return new{eltype(ics.ics)}(S, S0)
    end
end

# condition 
function (model::ImmortalModel)(u, t, integrator)
    return false
end

# affect! 
function (model::ImmortalModel)(integrator)
    return nothing
end

"""
    struct BrooksModelParameters{I, U, T, F}

A container for the parameters of the model of [Brooks et al. (2018)](https://www.int-res.com/abstracts/meps/v599/p1-18/).

### Parameters 

- `temp` [°C]: An [`InterpolatedField`](@ref) for the water temperature. Default `TEMPERATURE_ITP.x`.
- `no3` [mmol N/m^3]: An [`InterpolatedField`](@ref) for the Nitrogen content of the water. `NUTRIENTS_ITP.x`.
- `μ_max` [1/d]: Sargassum maximum growth rate. Value: `0.1`
- `m` [1/d]: Sargassum mortality rate. Value: `0.05`
- `I_k` [W/m^2]: Sargassum growth-irradiance parameter. Value: `70.0`
- `a_ref` [d]: Reference age for Sargassum light limitation. Value: `55.0`
- `k_N` [mmol N/m^3]: Sargassum nutrient (N) uptake half saturation. Value: `0.012`
- `T_ref` [°C]: Minimum temperature for Sargassum growth. Value: `18.0`
- `z_max` [m]: Maximum depth before Sargassum buoyancy is compromised. Value: `120.0`.
- `clumps_limits`: A `Tuple` of the form `(n_clumps_min, n_clumps_max)`. These impose hard lower \
and upper limits on the total number of clumps that can exist at any specific time (the total number \
of clumps that can have ever existed - i.e. `n_clumps_tot` of [`RaftParameters`](@ref) - may be higher.) Default: `(0, 3000)`.
- `dSdt`: Compute the rate of change of the "amount" `S` according to the Brooks model.

### dSdt

This function is of the form `dSdt = growth_factors  - death_factors`.
- `growth_factors = μ_max * light_factor * age_factor * temperature_factor * nutrients_factor`
- `death_factors = m`

### Constructors 

The function `BrooksModelParameters(; parameters...)` is provided, where each parameters is a kwarg 
with the default values given above.
"""
struct BrooksModelParameters{I<:InterpolatedField, U<:Integer, T<:Real, F<:Function} 
    temp::I
    no3::I
    μ_max::T
    m::T
    I_k::T
    a_ref::T
    k_N::T
    T_ref::T
    z_max::T
    clumps_limits::Tuple{U, U}
    dSdt::F

    function BrooksModelParameters(;
        temp::I = TEMPERATURE_ITP.x,
        no3::I = NUTRIENTS_ITP.x,
        μ_max::Real = 0.1, 
        m::Real = 0.05,
        I_k::Real = 70.0,
        a_ref::Real = 55.0,
        k_N::Real = 0.012,
        T_ref::Real = 18.0,
        z_max::Real = 120.0,
        clumps_limits::Tuple{Integer, Integer} = (0, 3000)) where {I<:InterpolatedField}

        μ_max, m, I_k, a_ref, k_N, T_ref, z_max = promote(μ_max, m, I_k, a_ref, k_N, T_ref, z_max)

        function brooks_dSdt_clump(x::Real, y::Real, t::Real)
            light_factor = 1.0 # 1 - exp(I/I_k)
            age_factor = 1.0 # exp(-t/params.a_ref)
            temp_factor = temp.fields[:temp](x, y, t) > T_ref ? 1.0 : 0.0
            N_factor = 1.0/(k_N/no3.fields[:no3](x, y, t) + 1.0)
            return μ_max*light_factor*age_factor*temp_factor*N_factor - m
        end
    
        return new{
            typeof(temp), 
            eltype(clumps_limits), 
            typeof(μ_max), 
            typeof(brooks_dSdt_clump)}(temp, no3, μ_max, m, I_k, a_ref, k_N, T_ref, z_max, clumps_limits, brooks_dSdt_clump)
    end
end


"""
    mutable struct BrooksModel{B, U, T}

The growth/death model of [Brooks et al. (2018)](https://www.int-res.com/abstracts/meps/v599/p1-18/). 

### Fields 

- `params`: The [`BrooksModelParameters`](@ref) parameters of the model.
- `growths`:A `Vector` of indices of clumps that are to be grown (if any).
- `deaths`: A `Vector` of indices of clumps that are to be killed (if any).
- `verbose`: A `Bool` such that `verbose = true` will log times and labels of clumps that grow and die.

### Constructors

Use `BrooksModel(ics::InitialConditions; S0 = [1.0], params = BrooksModelParameters(), verbose = false)`.

### Callbacks

Use [`cb_growth_death`](@ref) to create a `DiscreteCallback` suitable for use with `OrdinaryDiffEq.solve`. 

At each time step ...
"""
mutable struct BrooksModel{B<:BrooksModelParameters, U<:Integer, T<:Real} <: AbstractGrowthDeathModel
    S::Vector{T}
    S0::Vector{T}
    params::B
    growths::Vector{U}
    deaths::Vector{U}
    verbose::Bool

    function BrooksModel(
        ics::InitialConditions; 
        S0::Vector{T} = [1.0], 
        params::B = BrooksModelParameters(), 
        verbose = false) where {B<:BrooksModelParameters, T<:Real}

        S = fill(rand(S0), length(ics.ics))

        return new{B, Int64, T}(S, S0, params, Int64[], Int64[], verbose)
    end
end
    
# condition 
function (model::BrooksModel)(u, t, integrator)
    model.growths = Int64[]
    model.deaths = Int64[]

    for i = 1:n_clumps(u)
        rhs = model.params.dSdt(clump_i(u, i)..., t)*model.S[i]*(t - integrator.tprev)#*n_clumps(u)
        model.S[i] += rhs

        if model.S[i] < 0 && n_clumps(u) - length(model.deaths) >= model.params.clumps_limits[1]
            push!(model.deaths, i)
        elseif model.S[i] > 2 && n_clumps(u) + length(model.growths) <= model.params.clumps_limits[2]
            push!(model.growths, i)
        end
    end

    if length(model.growths) > 0 || length(model.deaths) > 0
        return true
    else
        return false
    end
end

# affect!
function (model::BrooksModel)(integrator)
    if model.verbose
        growths = [integrator.p.n_clumps_tot + i for i = 1:length(model.growths)]
        deaths = [integrator.p.loc2label[integrator.t][i] for i in model.deaths]
        if length(growths) > 0
            @info "Growth [t = $(integrator.t)]: $(growths)"
        end

        if length(deaths) > 0
            @info "Death [t = $(integrator.t)]: $(deaths)"
        end
    end

    ### growth/death specific
    for i in model.growths
        grow!(integrator, location = i)
    end

    kill!(integrator, model.deaths)

    ### random
    # for i in model.growths
    #     grow!(integrator, location = "parent")
    # end

    # for i = 1:length(model.deaths)
    #     kill!(integrator, rand(1:n_clumps(integrator.u)))
    # end

    return nothing
end