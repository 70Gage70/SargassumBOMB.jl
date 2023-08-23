using JLD2

include(joinpath(@__DIR__, "..", "interpolants", "interpolant-core.jl"))

########################################################################

# loading land interpolant and default EquirectangularReference
itp_path = joinpath(@__DIR__, "..", "interpolants", "land")
isdefined(@__MODULE__, :land_itp) || (const land_itp = load(joinpath(itp_path, "land_itp.jld2"), "land_itp"))
isdefined(@__MODULE__, :ref_itp) || (const ref_itp = land_itp.ref)

"""
    mutable struct Land{I, U}

A container for data handling death of clumps upon reaching the shore.

### Fields 

- `land_itp`: A `InterpolatedField` such that `land_itp.fields[:land](x, y)` is equal to `1.0` if `(x, y)` is on land and `0.0` otherwise.
- `deaths`: A `Vector` of indices of clumps that are to be killed.

### Constructors 

Use `Land(land_itp::InterpolatedField)` to create a new `Land` object.

### Callbacks 

Use `callback(land::Land)` to create a `DiscreteCallback` suitable for use with `OrdinaryDiffEq.solve`. At each time step, the position of 
each clump is checked and clumps that are currently on land are killed with [`kill!`](@ref).
"""
mutable struct Land{I<:InterpolatedField, U<:Integer}
    land_itp::I
    deaths::Vector{U}

    function Land(land_itp::InterpolatedField)
        return new{typeof(land_itp), Int64}(land_itp, Int64[])
    end
end

# condition 
function (land::Land)(u, t, integrator)
    land.deaths = [i for i = 1:n_clumps(u) if land.land_itp.fields[:land](clump_i(u, i)...) == 1.0]
    return length(land.deaths) > 0
end

# affect!
function (land::Land)(integrator)
    kill!(integrator, land.deaths)
end

# callback 
function callback(land::Land)
    return DiscreteCallback(land, land)
end