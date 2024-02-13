"""
    ymw2time(y, m, w)

Convert the time corresponding to the year `y`, month `m` and week `w` indicated into a single time measured in
days since `WATER_ITP.x.time_start`.

The days of the four weeks per month are defined as the 7th, 14th, 21nd and 28th.

### Example
```julia-repl
julia> WATER_ITP.x.time_start
2018-01-01T00:00:00

julia> ymw2time(2018, 1, 2)
14.0
```

Note that `t = 14.0` is midnight on January 15th, i.e. integrating from `t = 0.0` to `t = 14.0` will include the events of January 14th.
"""
function ymw2time(year::Integer, month::Integer, week::Integer)

    @assert week in [1, 2, 3, 4]

    return 1 + Day(DateTime(year, month, 7*week) - WATER_ITP.x.time_start).value |> float
end

"""
    ymwspan2weekspan(ymw1, ymw2)

Convert a time span between `ym1 = (year1, month1, week1)` and `ym2 = (year2, month2, week2)` to
a vector giving the time slices (measured in days since `WATER_ITP.x.time_start`) of each week.

The days of the four weeks per month are defined as the 7th, 14th, 21nd and 28th.

### Example
"""
function ymwspan2weekspan(ymw1::NTuple{3, Integer}, ymw2::NTuple{3, Integer})

    @assert ymw1[3] in [1, 2, 3, 4]
    @assert ymw2[3] in [1, 2, 3, 4]
    @assert ymw1 < ymw2

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


    return [(ymw2time(times[i]...), ymw2time(times[i + 1]...)) for i = 1:length(times) - 1]
end