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