"""
    T_REF

The time to which all times are referred.

This is a `Ref`, access or modify the actual value with `T_REF.x`.
"""
const T_REF = Ref{DateTime}(DateTime(2000, 1, 1))

"""
    datetime2time(dt)

Convert `dt::DateTime` to the amount of time since [`T_REF`](@ref) expressed in the units of `UNITS["time"]`.

This inverts [`time2datetime`](@ref).
"""
function datetime2time(dt::DateTime)
    return uconvert(UNITS["time"], dt - T_REF.x) |> x -> float(x.val)
end

"""
    time2datetime(time)

Convert the amount of time since [`T_REF`](@ref) expressed in the units of `UNITS["time"]` to a `DateTime`.

By convention, `time` is rounded down to the nearest second.   

This inverts [`datetime2time`](@ref).
"""
function time2datetime(time::Real)
    t = uconvert(u"s", time * UNITS["time"]) |> x -> round(Int64, x.val)

    return T_REF.x + Second(t)
end

"""
    ymw2time(y, m, w)

Convert the time corresponding to the year `y`, month `m` and week `w` indicated into a single time measured in
days since [`T_REF`](@ref).

The days of the four weeks per month are defined as the 7th, 14th, 21nd and 28th.
"""
function ymw2time(year::Integer, month::Integer, week::Integer)

    @argcheck week in [1, 2, 3, 4]

    return datetime2time(DateTime(year, month, 7*week))
end

"""
    ymwspan2weekspan(ymw1, ymw2)

Return a vector list of all `(year, month, week)` tuples between `ym1 = (year1, month1, week1)` and `ym2 = (year2, month2, week2)` inclusive.

### Example
"""
function ymwspan2weekspan(ymw1::NTuple{3, Integer}, ymw2::NTuple{3, Integer})

    @argcheck ymw1[3] in [1, 2, 3, 4]
    @argcheck ymw2[3] in [1, 2, 3, 4]
    @argcheck ymw1 < ymw2

    times = [ymw1]
    cur_time = ymw1

    while cur_time < ymw2
        nexty, nextm, nextw = cur_time .+ (0, 0, 1)

        if nextw == 5
            nextw = 1
            nextm = nextm + 1

            if nextm == 13
                nextm = 1
                nexty = nexty + 1
            end
        end

        cur_time = (nexty, nextm, nextw)
        push!(times, cur_time)
    end


    return times
end