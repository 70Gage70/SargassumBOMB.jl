using Dates

#########################

"""
    reduce_vector_to_range(vector)

Turn a vector `vector` of linearly spaced values [x0, x0 + δ, x0 + 2δ, ..., x0 + Nδ] into an `StepRangeLen` of the form `x0:δ:x0 + Nδ`.

δ may be positive or negative, but not zero.

Note that `reduce_vector_to_range(collect(x0:δ:x0 + Nδ)) = x0:δ:x0 + Nδ`.
"""
function reduce_vector_to_range(vector::Vector{<:Real})
    delta = vector[2] - vector[1]

    @assert delta != 0.0 "The first two points on the grid are at the same location."

    for i = 1:length(vector) - 1
        if vector[i + 1] - vector[i] != delta
            error("The grid is not uniformly spaced at entry $i.")
        end
    end

    start, stop = (vector[1], vector[end])
    return start:delta:stop
end

"""
    rata2datetime_minute(days)

Convert the number of Rata Die `days` into a `DateTime` rounded to the nearest minute.
"""
function rata2datetime_minute(days::Real)
    partial_days, full_days = modf(days)
    partial_days = round(partial_days*24*60)
    return rata2datetime(full_days) + Minute(partial_days)
end