"""
    abstract type AbstractGrowthDeathModel

The abstract type for growth and death models.

Subtypes must have fields `growths::Vector{<:Integer}`, `deaths::Vector{<:Integer}` and a 
function `dSdt`. This function must be callable at the solution vector and current time, 
i.e. it must be a function `dSdt(u, t)` which returns a `Real`.

The first entry in the solution vector is integrated in the sense that `du[1]/dt = dSdt`; 
this sets the "target" for the number of clumps that are alive or dead.
"""
abstract type AbstractGrowthDeathModel end 

"""
    mutable struct ImmortalModel{U, F}

An `AbstractGrowthDeathModel` such that no growth or death occurs.

### Constructors 

Use `ImmortalModel().`
"""
mutable struct ImmortalModel{U<:Integer, F<:Function} <: AbstractGrowthDeathModel
    growths::Vector{U}
    deaths::Vector{U}
    dSdt::F

    function ImmortalModel()
        return new{Int64, Function}(Int64[], Int64[], (u, t) -> 0.0)
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
    struct BrooksModelParameters{I, U, T}

A container for the interpolants and parameters of the model of [Brooks et al. (2018)](https://www.int-res.com/abstracts/meps/v599/p1-18/).

### Interpolants 

- `temp` [°C]: An [`InterpolatedField`](@ref) for the water temperature.
- `no3` [mmol N/m^3]: An [`InterpolatedField`](@ref) for the NO3 content of the water.

### Parameters 

- `μ_max` [1/d]: Sargassum maximum growth rate. Value: `0.1`
- `m` [1/d]: Sargassum mortality rate. Value: `0.05`
- `I_k` [W/m^2]: Sargassum growth-irradiance parameter. Value: `70.0`
- `a_ref` [d]: Reference age for Sargassum light limitation. Value: `55.0`
- `k_N` [mmol N/m^3]: Sargassum nutrient (N) uptake half saturation. Value: `0.012`
- `T_ref` [°C]: Minimum temperature for Sargassum growth. Value: `18.0`
- `z_max` [m]: Maximum depth before Sargassum buoyancy is compromised. Value: `120.0`.
- `clumps_limits`: A `Tuple` of the form `(n_clumps_min, n_clumps_max)`. These impose hard lower \
and upper limits on the total number of clumps that can exist at any specific time (the total number \
of clumps that can have ever existed - i.e. `n_clumps_tot` of [`RaftParameters`](@ref) - may be higher.) Default: `(0, 10000)`.

### Constructors 

The function `BrooksModelParameters(temp, no3; parameters...)` is provided, where each parameters is a kwarg 
with the default values given above.
"""
struct BrooksModelParameters{I<:InterpolatedField, U<:Integer, T<:Real} 
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

    function BrooksModelParameters(
        temp::InterpolatedField,
        no3::InterpolatedField;
        μ_max::Real = 0.1, 
        m::Real = 0.05,
        I_k::Real = 70.0,
        a_ref::Real = 55.0,
        k_N::Real = 0.012,
        T_ref::Real = 18.0,
        z_max::Real = 120.0,
        clumps_limits::Tuple{Integer, Integer} = (0, 10000))

        μ_max, m, I_k, a_ref, k_N, T_ref, z_max = promote(μ_max, m, I_k, a_ref, k_N, T_ref, z_max)
    
        return new{typeof(temp), eltype(clumps_limits), typeof(μ_max)}(temp, no3, μ_max, m, I_k, a_ref, k_N, T_ref, z_max, clumps_limits)
    end
end

"""
    brooks_dSdt_clump(x, y, t; params::BrooksModelParameters, n_clumps)

Compute `dS/dt` for a single clump at position `(x, y)` and time `t` using the [`BrooksModelParameters`](@ref) in `params`. 

### Optional Arguments

- `n_clumps`: An `Integer` number of total clumps at the current time. Not present in the original Brooks model.
"""
function brooks_dSdt_clump(x::Real, y::Real, t::Real, params::BrooksModelParameters; n_clumps::Integer)
    light_factor = 1.0 # 1 - exp(I/I_k)
    age_factor = 1.0 # exp(-t/params.a_ref)
    temp_factor = params.temp.fields[:temp](x, y, t) > params.T_ref ? 1.0 : 0.0
    N_factor = 1.0/(params.k_N/params.no3.fields[:no3](x, y, t) + 1.0)
    return params.μ_max*light_factor*age_factor*temp_factor*N_factor - params.m
end

"""
    brooks_dSdt_raft(u, t, params::BrooksModelParameters)

Compute `dS/dt` for a raft with solution vector `u` at time `t`.

Calculates the median of [`brooks_dSdt_clump`](@ref) evaluated at each clump multiplied by `u[1]`, the "amount" parameter in the solution vector.
"""
function brooks_dSdt_raft(u::Vector{<:Real}, t::Real, params::BrooksModelParameters)
    return median(brooks_dSdt_clump(clump_i(u, i)..., t, params, n_clumps = n_clumps(u)) for i = 1:n_clumps(u))*u[1]
end

"""
    mutable struct BrooksModel{B, U, F}

The growth/death model of [Brooks et al. (2018)](https://www.int-res.com/abstracts/meps/v599/p1-18/). 

### Fields 

- `params`: The [`BrooksModelParameters`](@ref) parameters of the model.
- `dSdt`: A `Function` mapping `(u, t)` to `dS/dt`, given by [`brooks_dSdt_raft`](@ref).
- `growths`:A `Vector` of indices of clumps that are to be grown.
- `deaths`: A `Vector` of indices of clumps that are to be killed.
- `verbose`: A `Bool` such that `verbose = true` will log times and labels of clumps that grow and die.

### Constructors

Use `BrooksModel(;params = BrooksModelParameters(TEMPERATURE_ITP, NUTRIENTS_ITP), verbose = false)`.

### Callbacks

Use [`cb_growth_death`](@ref) to create a `DiscreteCallback` suitable for use with `OrdinaryDiffEq.solve`. At each time step, the 
value of `dSdt` is evaluated and clumps are grown with [`grow!`](@ref) and killed with [`kill!`](@ref) according to the following 
default logic:

- When the difference between `u[1]` and the actual number of clumps is at least 1 (call it `δn`), the `δn` clumps with the most \
extreme values of [`brooks_dSdt_clump`](@ref) are selected. If `δn < 0`, they are killed. If `δn > 0`, then those clumps are chosen \
as parents in the [`grow!`](@ref) method.
"""
mutable struct BrooksModel{B<:BrooksModelParameters, U<:Integer, F<:Function} <: AbstractGrowthDeathModel
    params::B
    dSdt::F
    growths::Vector{U}
    deaths::Vector{U}
    verbose::Bool

    function BrooksModel(;params::BrooksModelParameters = BrooksModelParameters(TEMPERATURE_ITP.x, NUTRIENTS_ITP.x), verbose = false)
        return new{typeof(params), Int64, Function}(params, (u, t) -> brooks_dSdt_raft(u, t, params), Int64[], Int64[], verbose)
    end
end
    
# condition 
function (model::BrooksModel)(u, t, integrator)
    delta_n = round(Int64, u[1]) - n_clumps(u)
    delta_n = sign(delta_n)*min(round(Int64, 0.05*n_clumps(u)), abs(delta_n))
    # ensures that delta_n does not change the number of clumps by more than 5% at any step

    if (delta_n > 0 && u[1] + delta_n > model.params.clumps_limits[2])
        delta_n = max(0, model.params.clumps_limits[2] - round(Int64, u[1]))
    elseif (delta_n < 0 && u[1] + delta_n < model.params.clumps_limits[1])
        delta_n = min(0, model.params.clumps_limits[1] - round(Int64, u[1]))
    end

    if delta_n != 0
        dSdt = [brooks_dSdt_clump(clump_i(u, i)..., t, model.params, n_clumps = n_clumps(u)) for i = 1:n_clumps(u)]

        if delta_n > 0 
            # need to grow clumps, take the delta_n clumps with largest positive dSdt as parents
            model.deaths = typeof(model.deaths)[]
            model.growths = partialsortperm(dSdt, 1:delta_n, rev = true) |> collect
        else 
            # need to kill clumps, take the delta_n with largest negative dSdt to kill
            model.deaths = partialsortperm(dSdt, 1:min(-delta_n, n_clumps(u))) |> collect
            model.growths = typeof(model.growths)[]
        end

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

    for i in model.growths
        # grow!(integrator, location = i)
        grow!(integrator, location = "parent")
    end

    # kill!(integrator, model.deaths)
    for i = 1:length(model.deaths)
        kill!(integrator, rand(1:n_clumps(integrator.u)))
    end

    return nothing
end