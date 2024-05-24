"""
    struct InitialConditions

A container for the initial conditions for a raft. 

### Fields

- `tspan`: A `Tuple` such that the integration is performed for `tspan[1] ≤ t ≤ tspan[2]` where `t` is \
in `UNITS["time"]` since [`T_REF`](@ref).
- `ics`: A `2 x N` `Matrix` of the form `[x1 x2 ... xN ; y1 y2 ... yN]` giving the initial coordinates of each clump.

### Constructors 

use `InitialConditions(;tspan, ics)`.
"""
struct InitialConditions
    tspan::Tuple{Float64, Float64}
    ics::Matrix{Float64}

    function InitialConditions(;tspan::Tuple{Real, Real}, ics::Matrix{<:Real})
        @argcheck tspan[1] < tspan[2] "initial time must be less than final time"
        @argcheck size(ics, 1) == 2 "matrix should have size 2 x N"
        @argcheck size(ics, 2) > 0 "ics can not be empty"

        return new(tspan, ics)
    end
end

"""
    InitialConditions(tspan, xy0; to_xy)

Construct initial conditions suitable for use in `RaftParameters.ics` from `2 x N` `Matrix`, `xy0` which should \
be equirectangular coordinates

Can be applied as `InitialConditions(tspan, x_range, y_range; [kwargs])` to generate clumps in a rectangular arrangement.

Can be applied as `InitialConditions(tspan, x0, y0; to_xy; [kwargs])` for a single clump with coordinates `(x0, y0)`.

### Optional Arguments 

`to_xy`: If `true`, the coordinates are converted from spherical to equirectangular coordinates. Default `false`.
"""
function InitialConditions(
    tspan::Tuple{Real, Real}, 
    xy0::Matrix{<:Real};
    to_xy::Bool = false)
    
    if to_xy
        ics = sph2xy(xy0)
    else
        ics = deepcopy(xy0)
    end

    return InitialConditions(tspan = tspan, ics = ics)
end

function InitialConditions(
    tspan::Tuple{Real, Real},
    x_range::AbstractRange{T}, 
    y_range::AbstractRange{T};
    to_xy::Bool = false) where {T<:Real}

    @argcheck allunique(x_range) "`x_range` can not have repeated entries"
    @argcheck allunique(y_range) "`y_range` can not have repeated entries"

    if to_xy
        ics_x, ics_y = sph2xy(x_range, y_range)
    else
        ics_x, ics_y = x_range, y_range
    end

    ics = Iterators.product(ics_x, ics_y) |> collect |> vec |> stack

    return InitialConditions(tspan = tspan, ics = ics)
end

function InitialConditions(
    tspan::Tuple{Real, Real}, 
    x0::Real, 
    y0::Real;
    to_xy::Bool = false)

    if to_xy
        ics = sph2xy(x0, y0)
    else
        ics = [x0, y0]
    end

    return InitialConditions(tspan = tspan, ics = reshape(ics, 2, 1))
end

"""
    InitialConditions(tspan, dist, weeks, number, sample_type; n_clumps_max, seed)

Construct [`InitialConditions`](@ref) from a `SargassumDistribution`.

### Arguments 
- `tspan`: A `Tuple{Real, Real}` such that the integration is performed for `tspan[1] ≤ t ≤ tspan[2]` where `t` is \
in in `UNITS["time"]` since [`T_REF`](@ref).
- `dist`: A `SargassumDistribution`.
- `weeks`: A `Vector{<:Integer}` giving the weeks of the month to consider. Each entry should be between 1 and 4 and appear only once.
- `number`: The number of clumps to initialize; interactive with `sample_type` and should be at least `1`.
- `sample_type`: A `String` identifying one of the methods of assigning clump locations based on the distribution.
    - `"levels"`: Boxes with nonzero Sargassum are divided into `number` levels of size `(minimum(dist.sargassum) - maximum(dist.sargassum))/number`.
    Each box gets a number of clumps equal to its level index. For example, if `number = 2`, then the smaller half of the boxes (by Sargassum content) 
    get 1 clump each and the larger half get 2 clumps each.
    - `"sample"`: A number `number` of samples are drawn from `dist`. Each sample is placed uniformly at random inside the corresponding box.
    - `"sorted"`: Boxes are filled with one clump placed uniformly at random inside them, starting from the box with the highest concentration. If `number` 
                is greater than the total number of boxes, repeat the loop starting again from the highest concentration box.
    - `"uniform"`: Exactly one clump is placed in the center of each box with nonzero concentration. Note that this ignores `number`.

### Optional Arguments 

- `seed`: A `Random.seed!` used in the initialization.
"""
function InitialConditions(
    tspan::Tuple{Real, Real},
    dist::SargassumDistribution, 
    weeks::Vector{<:Integer},
    number::Integer,
    sample_type::String;
    seed::Integer = 1234)

    seed!(seed)

    @argcheck sample_type in ["sample", "sorted", "uniform", "levels"] "`sample_type` not recognized."
    @argcheck number > 0 "Must request at least one clump"
    @argcheck length(weeks) > 0 "At least one week must be selected."
    @argcheck all(map(x -> 1 <= x <= 4, weeks)) "Each entry of `weeks` must be between 1 and 4."
    @argcheck allunique(weeks) "Each week should only appear once."

    sarg = sum(dist.sargassum[:,:,week] for week in weeks)
    lons = dist.lon
    lats = dist.lat

    δ_x = abs(lons[2] - lons[1])
    δ_y = abs(lats[2] - lats[1])

    x0 = Float64[]
    y0 = Float64[]

    pts = Iterators.product(lons, lats) |> x -> reduce(vcat, collect(x)) # vector of (lon, lat)
    wts = sarg |> x -> reduce(vcat, x) # vector of sarg, matched with pts

    if sample_type == "sample"
        samples = Distributions.sample(pts, Weights(wts), number)

        for samp in samples
            push!(x0, rand(Uniform(samp[1] - δ_x/2, samp[1] + δ_x/2)))
            push!(y0, rand(Uniform(samp[2] - δ_y/2, samp[2] + δ_y/2)))
        end
    elseif sample_type == "sorted"
        idx = findall(x -> x > 0, wts)
        pts = pts[idx]
        wts = wts[idx]

        pts = pts[sortperm(wts, rev = true)]

        n_c = 1
        while n_c <= number
            idx = mod(n_c, length(pts)) == 0 ? length(pts) : mod(n_c, length(pts))
            pt = pts[idx]

            push!(x0, rand(Uniform(pt[1] - δ_x/2, pt[1] + δ_x/2)))
            push!(y0, rand(Uniform(pt[2] - δ_y/2, pt[2] + δ_y/2)))       

            n_c = n_c + 1
        end
    elseif sample_type == "uniform"
        idx = findall(x -> x > 0, wts)
        pts = pts[idx]

        for pt in pts
            push!(x0, pt[1])
            push!(y0, pt[2])  
        end
    elseif sample_type == "levels"
        idx = findall(x -> x > 0, wts)
        pts = pts[idx]
        wts = wts[idx]

        wtsmin, wtsmax = extrema(wts)

        for i = 1:length(pts)
            pt = pts[i]
            n_c = number*(wts[i] - wtsmin)/(wtsmax - wtsmin) |> x -> round(Integer, x, RoundUp)
            for _ = 1:n_c
                push!(x0, rand(Uniform(pt[1] - δ_x/2, pt[1] + δ_x/2)))
                push!(y0, rand(Uniform(pt[2] - δ_y/2, pt[2] + δ_y/2)))       
            end
        end           
    end

    ics = [x0' ; y0']
    ics = sph2xy(ics)
    
    return InitialConditions(tspan = tspan, ics = ics)
end