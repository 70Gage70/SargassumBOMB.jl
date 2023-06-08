using MAT

##############################################

abstract type AbstractVectorField end


mutable struct VectorField2DGrid{T<:Real} <: AbstractVectorField
    lon::AbstractRange{T}
    lat::AbstractRange{T}
    time0::T
    time::AbstractRange{T}
    u::Array{T,3}
    v::Array{T,3}
end


function reduce_to_range(vector::Vector{<:Real})
    delta = vector[2] - vector[1]

    @assert delta != 0.0 "The first two points on the grid are at the same location."

    for i = 1:length(vector) - 1
        if vector[i + 1] - vector[i] != delta
            error("The grid is not uniformly spaced at entry $i.")
        end
    end

    start, stop = extrema(vector)
    return start:delta:stop
end


function VectorField2DGrid(
    infile::String;
    lon_alias::String = "lon",
    lat_alias::String = "lat",
    time_alias::String = "t",
    u_alias::String = "u",
    v_alias::String = "v",
    NaN_replacement::Any = 0.0)

    # loading data
    extension = infile[findlast(==('.'), infile)+1:end]
    @assert extension in ["mat"] "Require a .mat file."

    data = matopen(infile)

    @assert lon_alias in collect(keys(data)) "$(lon_alias) not in $(infile)"
    @assert lat_alias in collect(keys(data)) "$(lat_alias) not in $(infile)"
    @assert time_alias in collect(keys(data)) "$(time_alias) not in $(infile)"
    @assert u_alias in collect(keys(data)) "$(u_alias) not in $(infile)"
    @assert v_alias in collect(keys(data)) "$(v_alias) not in $(infile)"

    lon, lat, time = map(vec, read(data, lon_alias, lat_alias, time_alias));
    u, v = read(data, u_alias, v_alias);
    close(data)

    # matching dimensions of lon, lat time with u and v
    # also ensure that u is of the form u[lon, lat, time] and the same for v
    var_sizes = [length(var) for var in [lon, lat, time]]
    u_size = size(u)
    v_size = size(v)
    
    perm_u = Vector{Int64}()
    for size in var_sizes
        dim = findfirst(x->x==size, u_size)
        @assert dim!==nothing "The dimensions of the variables and u vector field components must match. Got $(size) and $(u_size)."

        push!(perm_u, dim)
    end

    perm_v = Vector{Int64}()
    for size in var_sizes
        dim = findfirst(x->x==size, v_size)
        @assert dim!==nothing "The dimensions of the variables and v vector field components must match. Got $(size) and $(v_size)."

        push!(perm_v, dim)
    end    

    u = permutedims(u, perm_u)
    v = permutedims(v, perm_v)

    # ensure that lon, lat, time are defined on an (increasing) regularly spaced grid
    lon = reduce_to_range(lon)
    lat = reduce_to_range(lat)
    time = reduce_to_range(time)

    if step(lon) < 0
        reverse!(u, dims = 1)
        reverse!(v, dims = 1)
        lon = reverse(lon)
    end

    if step(lat) < 0
        reverse!(u, dims = 2)
        reverse!(v, dims = 2)
        lat = reverse(lat)
    end

    if step(time) < 0
        reverse!(u, dims = 3)
        reverse!(v, dims = 3)
        time = reverse(time)
    end

    # NaN replacement
    u[isnan.(u)] .= NaN_replacement
    v[isnan.(v)] .= NaN_replacement

    # rescale time
    time0 = first(time)
    time = time .- time0

    return VectorField2DGrid(lon, lat, time0, time, u, v)
end