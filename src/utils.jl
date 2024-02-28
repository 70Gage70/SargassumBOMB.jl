"""
    vec2range(vector)

Convert a `Vector` of linearly spaced values to a `StepRangeLen`.
"""
function vec2range(vector::Vector{<:Real})
    step_size = vector[2] - vector[1]
    for i in 2:length(vector)-1
        if abs(vector[i+1] - vector[i] - step_size) > 1e-10
            error("Input vector is not linearly spaced")
        end
    end
    
    return vector[1]:step_size:vector[end]
end