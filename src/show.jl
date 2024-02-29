function Base.show(io::IO, x::EquirectangularReference)
    print(io, "EquirectangularReference[lon0 = $(x.lon0), lat0 = $(x.lat0), R = $(x.R)]")
end

function Base.show(io::IO, x::ClumpParameters)
    α, τ, R, f, σ = (x.α, x.τ, x.R, x.f, x.σ) .|> x-> round(x, sigdigits = 4)
    print(io, "ClumpParameters[α = $α, τ = $τ, R = $R, f = $f, σ = $σ]")
end

function Base.show(io::IO, x::SpringParameters)
    print(io, "SpringParameters[1->")
    show(io, x.k(1))
    print(io, ", ")
    show(io, x.L)
    print(io, "]")
end

function Base.show(io::IO, x::InitialConditions)
    print(io, "InitialConditions[")
    print(io, "time ∈ ($(time2datetime(x.tspan[1])), $(time2datetime(x.tspan[2]))), ")
    print(io, "n_clumps = $(x.ics[1]), ")
    xmin, xmax = extrema(x.ics[2:2:end])
    ymin, ymax = extrema(x.ics[3:2:end])
    lon_min, lat_min = xy2sph(xmin, ymin) .|> x -> round(x, sigdigits = 4)
    lon_max, lat_max = xy2sph(xmax, ymax) .|> x -> round(x, sigdigits = 4)
    print(io, "lon/lat ∈ ($lon_min, $lon_max) × ($lat_min, $lat_max)]")
end

function Base.show(io::IO, x::RaftParameters)
    println(io, " RaftParameters")
    println(io, "ICS = $(x.ics)")
    println(io, "Clumps = $(x.clumps)")
    println(io, "Springs = $(x.springs)")
    println(io, "Connections = $(typeof(ConnectionsNearest(3)).name.name)")
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
    println(io, " temp = TEMPERATURE_ITP.x")
    println(io, " no3 = NUTRIENTS_ITP.x")
    println(io, " μ_max = $(x.μ_max)")
    println(io, " m = $(x.m)")
    println(io, " I_k = $(x.I_k)")
    println(io, " a_ref = $(x.a_ref)")
    println(io, " k_N = $(x.k_N)")
    println(io, " T_ref = $(x.T_ref)")
    println(io, " z_max = $(x.z_max)")
    println(io, " clumps_limits = $(x.clumps_limits)")
end

function Base.show(io::IO, x::BrooksModel)
    print(io, "BrooksModel")
end

function Base.show(io::IO, bop::BOMBOptimizationProblem)
    optimizable = [bop.params[param].name for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]
    space = [(name, bop.params[name].default, (bop.params[name].bounds)) for name in optimizable]
    optimized_q = bop.opt === nothing ? false : true
    
    println(io, "BOMBOptimizationProblem")

    println(io, " Optimized?: $(optimized_q)")
    println(io, " Search space: $(space)")

    if optimized_q
        opts = Tuple{String, Float64}[]
        for name in OPTIMIZATION_PARAMETER_NAMES
            opt = bop.params[name].opt
            if opt === nothing 
                opt = bop.params[name].default
            end
            push!(opts, (name, round(Float64(opt), sigdigits = 4)))
        end
        println(io, " Optimals: $(opts)")
        println(io, " Loss: $(bop.opt)")
    end
end