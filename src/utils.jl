# """
#     n_clumps(object)

# Return the number of living clumps in the `object.`

# `object` can be a [`InitialConditions`](@ref). This simply returns `size(ics.ics, 2)`.

# `object` can be an [`RaftParameters`](@ref).
# """
# function n_clumps(ics::InitialConditions)
#     return size(ics.ics, 2)
# end

"""
    clump_i(u, i)

Return a view to the the `[x, y]` coordinates of the `i`th clump in the solution matrix `u`. This is `view(u :,i)`.
"""
function clump_i(u::Matrix{Float64}, i::Integer)
    return view(u, :, i)
end

"""
    com(u)

Return the center of mass `[x, y]` coordinates of the solution matrix `u`.
"""
function com(u::Matrix{Float64})
    return [mean(u[1,:]), mean(u[2,:])]
end

"""
    vec2range(vector; force)

Convert a `Vector` of linearly spaced values to a `StepRangeLen`. 

If `force == true`, the range will be constructed even if the vector isn't linearly spaced by linear interpolation preserving the length. Default `false`.
"""
function vec2range(vector::Vector{<:Real}; force::Bool = false)
    if force
        return range(vector[1], vector[end], length = length(vector))
    end

    step_size = vector[2] - vector[1]
    for i in 2:length(vector)-1
        if abs(vector[i+1] - vector[i] - step_size) > 1e-10
            error("Input vector is not linearly spaced, expected $(step_size), got $(abs(vector[i+1] - vector[i])) at positions $((i, i + 1))")
        end
    end
    
    return vector[1]:step_size:vector[end]
end