"""
    n_clumps(u)

Return the number of clumps in the solution vector `u`. This is `floor(Int64, length(u)/2)`.
"""
function n_clumps(u::Vector{<:Real})
    return floor(Int64, length(u)/2)
end

"""
    clump_i(u, i)

Return the `[x, y]` coordinates of the `i`th clump in the solution vector `u`. This is `u[2*i-1:2*i]`.
"""
function clump_i(u::Vector{<:Real}, i::Integer)
    return u[2*i-1:2*i]
end

"""
    com(u)

Return the center of mass `[x, y]` coordinates of the solution vector `u`.
"""
function com(u::Vector{<:Real})
    return [mean(u[1:2:end]), mean(u[2:2:end])]
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