function Base.show(io::IO, x::EquirectangularReference)
    print(io, "EquirectangularReference[lon0 = $(x.lon0), lat0 = $(x.lat0), R = $(x.R)]")
end

function Base.show(io::IO, x::Union{GriddedField, InterpolatedField})
    if x isa GriddedField
        println(io, "GriddedField")
    else
        println(io, "InterpolatedField")
    end

    print(io, "Dimensions = ")
    for y in x.dims_names
        print("$(y) ")
    end

    println()    

    print(io, "Fields = ")
    for y in x.fields_names
        print("$(y) ")
    end

    println()

    if (:x in keys(x.dims)) && (:y in keys(x.dims))
        lon, lat = xy2sph(x.dims[:x], x.dims[:y]) .|> extrema
        lon = lon .|> x -> round(x, sigdigits = 4)
        lat = lat .|> x -> round(x, sigdigits = 4)
        println(io, "lon/lat ∈ ($(lon[1]), $(lon[2])) × ($(lat[1]), $(lat[2]))]")
    end    

    if (:t in keys(x.dims))
        tmin, tmax = extrema(x.dims[:t])
        print(io, "time ∈ ($(time2datetime(tmin)), $(time2datetime(tmax)))")
    end
end

function Base.show(io::IO, x::ClumpParameters)
    α, τ, R, f, σ = (x.α, x.τ, x.R, x.f, x.σ) .|> x-> round(x, sigdigits = 4)
    print(io, "ClumpParameters[α = $α, τ = $τ, R = $R, f = $f, σ = $σ]")
end

function Base.show(io::IO, x::HookeSpring)
    print(io, "HookeSpring[k = ")
    show(io, x.k)
    print(io, ", L = ")
    show(io, x.L)
    print(io, "]")
end

function Base.show(io::IO, x::BOMBSpring)
    print(io, "BOMBSpring[A = ")
    show(io, x.A)
    print(io, ", L = ")
    show(io, x.L)
    print(io, "]")
end

function Base.show(io::IO, x::InitialConditions)
    print(io, "InitialConditions[")
    print(io, "time ∈ ($(time2datetime(x.tspan[1])), $(time2datetime(x.tspan[2]))), ")
    print(io, "n_clumps = $(size(x.ics, 2)), ")
    xmin, xmax = extrema(x.ics[1:2:end])
    ymin, ymax = extrema(x.ics[2:2:end])
    lon_min, lat_min = xy2sph(xmin, ymin) .|> x -> round(x, sigdigits = 4)
    lon_max, lat_max = xy2sph(xmax, ymax) .|> x -> round(x, sigdigits = 4)
    print(io, "lon/lat ∈ ($lon_min, $lon_max) × ($lat_min, $lat_max)]")
end

function Base.show(io::IO, x::RaftParameters)
    println(io, " RaftParameters")
    println(io, "ICS = $(x.ics)")
    println(io, "Clumps = $(x.clumps)")
    println(io, "Springs = $(x.springs)")
    println(io, "Connections = $(typeof(x.connections).name.name)")
    println(io, "GrowthDeath = $(x.gd_model)")
end

function Base.show(io::IO, x::Trajectory)
    if length(x.t) == 0
        print(io, "Trajectory[0 pts]")
    else
        print(io, "Trajectory[(")
        show(io, first(x.t))
        print(io, ", ")
        show(io, last(x.t))
        print(io, "), ")
        show(io, length(x.t))
        print(io, " pts]")
    end
end

function Base.show(io::IO, x::RaftTrajectory)
    print(io, "RaftTrajectory[")
    show(io, length(keys(x.trajectories)))
    print(io, " trajectories, ")
    show(io, length(x.t))
    print(io, " times]")
end

function Base.show(io::IO, x::NoLand)
    print(io, "NoLand")
end

function Base.show(io::IO, x::Land)
    print(io, "Land[land_itp = LAND_ITP.x]")
end

function Base.show(io::IO, x::ImmortalModel)
    print(io, "ImmortalModel")
end

function Base.show(io::IO, x::BrooksModelParameters)
    println(io, "BrooksModelParameters")
    println(io, " temp = $(x.temp)")
    println(io, " no3 = $(x.no3)")
    println(io, " μ_max = $(x.μ_max)")
    println(io, " m = $(x.m)")
    println(io, " I_k = $(x.I_k)")
    println(io, " a_ref = $(x.a_ref)")
    println(io, " k_N = $(x.k_N)")
    println(io, " T_ref = $(x.T_ref)")
    println(io, " z_max = $(x.z_max)")
    println(io, " clumps_limits = $(x.clumps_limits)")
    println(io, " dSdt = $(nameof(x.dSdt))")
end

function Base.show(io::IO, x::BrooksModel)
    print(io, "BrooksModel")
end

function Base.show(io::IO, x::TimeSeries)
    print(io, "TimeSeries[")
    print(io, "$(length(x.t)) weeks ∈ ($(x.t[1][1]), $(x.t[end][end])), ")
    lon_min, lon_max = extrema(x.lon) .|> x -> round(x, sigdigits = 4)
    lat_min, lat_max = extrema(x.lat) .|> x -> round(x, sigdigits = 4)
    print(io, "lon/lat ∈ ($lon_min, $lon_max) × ($lat_min, $lat_max), ") 
    if x.exclude_clouded_bins
        print(io, "clouds excluded]")
    else
        print(io, "clouds not excluded]")
    end
end