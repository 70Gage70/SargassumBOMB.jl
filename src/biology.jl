using OrdinaryDiffEq
using JLD2

include(joinpath(@__DIR__, "parameters.jl"))
include(joinpath(@__DIR__, "..", "interpolants", "itp-core.jl"))

########################################################################

# loading interpolants
itp_path = joinpath(@__DIR__, "..", "interpolants", "biology")
isdefined(@__MODULE__, :temp_itp) || (const temp_itp = load(joinpath(itp_path, "temp_itp.jld2"), "temp_itp"))
isdefined(@__MODULE__, :no3_itp) || (const no3_itp = load(joinpath(itp_path, "no3_itp.jld2"), "no3_itp"))

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

# callback 
function callback(model::ImmortalModel)
    return DiscreteCallback(model, model)
end

function Base.show(io::IO, x::ImmortalModel)
    print(io, "ImmortalModel")
end

"""
    struct BrooksModelParameters

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

### Constructors 

The function `BrooksModelParameters(temp, no3; constants...)` is provided, where each constant is a kwarg 
with the default values given above.
"""
struct BrooksModelParameters{I<:InterpolatedField, T<:Real} 
    temp::I
    no3::I
    μ_max::T
    m::T
    I_k::T
    a_ref::T
    k_N::T
    T_ref::T
    z_max::T


    function BrooksModelParameters(
        temp::InterpolatedField,
        no3::InterpolatedField;
        μ_max::Real = 0.1, 
        m::Real = 0.05,
        I_k::Real = 70.0,
        a_ref::Real = 55.0,
        k_N::Real = 0.012,
        T_ref::Real = 18.0,
        z_max::Real = 120.0)

        μ_max, m, I_k, a_ref, k_N, T_ref, z_max = promote(μ_max, m, I_k, a_ref, k_N, T_ref, z_max)
    
        return new{typeof(temp), typeof(μ_max)}(temp, no3, μ_max, m, I_k, a_ref, k_N, T_ref, z_max)
    end
end

"""
    dSdt1(x, y, t; params::BrooksModelParameters)
"""
function dSdt1(x::Real, y::Real, t::Real; params::BrooksModelParameters)
    light_factor = 1.0 # 1 - exp(I/I_k)
    age_factor = 1.0 # exp(-t/params.a_ref)
    temp_factor = params.temp.fields[:temp](x, y, t) > params.T_ref ? 1.0 : 0.0
    N_factor = 1.0/(params.k_N/params.no3.fields[:no3](x, y, t) + 1.0)
    return params.μ_max*light_factor*age_factor*temp_factor*N_factor - params.m
end

"""
    brooks_dSdt(u, t; params::BrooksModelParameters)
"""
function brooks_dSdt(u::Vector{<:Real}, t::Real; params::BrooksModelParameters)
    return mean(dSdt1(clump_i(u, i)..., t, params = params) for i = 1:n_clumps(u))*u[1]
end

"""
    mutable struct BrooksModel{B, T}
"""
mutable struct BrooksModel{B<:BrooksModelParameters, U<:Integer, F<:Function} <: AbstractGrowthDeathModel
    params::B
    growths::Vector{U}
    deaths::Vector{U}
    dSdt::F

    function BrooksModel(;params::BrooksModelParameters = BrooksModelParameters(temp_itp, no3_itp))
        return new{typeof(params), Int64, Function}(params, Int64[], Int64[], (u, t) -> brooks_dSdt(u, t, params = params))
    end
end
    
# condition 
function (model::BrooksModel)(u, t, integrator)
    delta_n = round(Int64, u[1]) - n_clumps(u)

    if delta_n != 0
        dSdt = [dSdt1(clump_i(u, i)..., t, params = model.params) for i = 1:n_clumps(u)]
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

    for i in model.growths
        # println("[BROOKS] Growing $([integrator.p.n_clumps_tot]) at time $(integrator.t)")
        grow!(integrator, locations = i, connections = "full")
    end

    # if length(model.deaths) > 0
    #     println("[BROOKS] Killing $(model.deaths) at time $(integrator.t)")
    # end

    kill!(integrator, model.deaths)

    return nothing
end

# callback 
function callback(brooks::BrooksModel)
    return DiscreteCallback(brooks, brooks)
end