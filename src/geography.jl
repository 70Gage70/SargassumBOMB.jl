"""
    abstract type AbstractLand

A supertype for all land/shore types.
"""
abstract type AbstractLand end 

"""
    struct NoLand

An [`AbstractLand`](@ref) such that the land/shore is completely ignored.
"""
struct NoLand <: AbstractLand
end

# condition 
function (land::NoLand)(u, t, integrator)
    return false
end

# affect!
function (land::NoLand)(integrator)
    return nothing
end

function Base.show(io::IO, x::NoLand)
    print(io, "NoLand")
end

"""
    mutable struct Land{I, U}

A container for data handling death of clumps upon reaching the shore.

### Fields 

- `land_itp`: A `InterpolatedField` such that `land_itp.fields[:land](x, y)` is equal to `1.0` if `(x, y)` is on land and `0.0` otherwise.
- `deaths`: A `Vector` of indices of clumps that are to be killed.
- `verbose`: A `Bool` such that `verbose = true` will log times and labels of clumps that hit land.

### Constructors 

Use `Land(;land_itp::InterpolatedField = land_itp, verbose = false)` to create a new `Land` object.
"""
mutable struct Land{I<:InterpolatedField, U<:Integer} <: AbstractLand
    land_itp::I
    deaths::Vector{U}
    verbose::Bool

    function Land(;land_itp::InterpolatedField = LAND_ITP.x, verbose = false)
        return new{typeof(land_itp), Int64}(land_itp, Int64[], verbose)
    end
end

# condition 
function (land::Land)(u, t, integrator)
    land.deaths = [i for i = 1:n_clumps(u) if land.land_itp.fields[:land](clump_i(u, i)...) == 1.0]
    return length(land.deaths) > 0
end

# affect!
function (land::Land)(integrator)
    if land.verbose
        deaths = [integrator.p.loc2label[integrator.t][i] for i in land.deaths]
        @info "Land [t = $(integrator.t)]: $(deaths)"
    end
    
    kill!(integrator, land.deaths)
end

function Base.show(io::IO, x::Land)
    print(io, "Land[land_itp = LAND_ITP.x]")
end