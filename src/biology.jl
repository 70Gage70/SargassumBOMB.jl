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
    brooks_dSdt(x, y, t; p::BrooksModelParameters)
"""
function brooks_dSdt(x::Real, y::Real, t::Real; params::BrooksModelParameters)
    temp_factor = params.temp.fields[:temp](x, y, t) > params.T_ref ? 1.0 : 0.0
    return params.μ_max*1.0*exp(-t/params.a_ref)*temp_factor*(1.0/(params.k_N/params.no3.fields[:no3](x, y, t) + 1.0)) - params.m
end

"""
    mutable struct BrooksModel{B, U}
"""
mutable struct BrooksModel{B<:BrooksModelParameters, U<:Integer}
    params::B
    deaths::Vector{U}
    growths::Vector{U}

    function BrooksModel(params::BrooksModelParameters)
        return new{typeof(params), Int64}(params, Int64[], Int64[])
    end
end
    
# condition 
function (model::BrooksModel)(u, t, integrator)
    dSdt = [brooks_dSdt(clump_i(u, i)..., t, params = model.params) for i = 1:n_clumps(u)]

    println("t = ", t, ": +", length(dSdt[dSdt .> 0.0]), " | -", length(dSdt[dSdt .< 0.0]))
    return false
end

# affect!
function (model::BrooksModel)(integrator)
    println("affecting!")
end

# callback 
function callback(brooks::BrooksModel)
    return DiscreteCallback(brooks, brooks)
end