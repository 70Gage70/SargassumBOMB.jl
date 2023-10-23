include(joinpath(@__DIR__, "SargassumBOMB.jl"))

using Optimization 

########################################################################

function integrate_water(
    α::Real, 
    τ::Real;
    initial_time::NTuple{2, Integer} = (2018, 3), 
    final_time::NTuple{2, Integer} = (2018, 4))

    tstart = Day(DateTime(initial_time...) - DateTime(yearmonth(water_itp.time_start)...)).value |> float
    tend = tstart + Day(DateTime(final_time...) - DateTime(initial_time...)).value
    tspan = (tstart, tend)

    dists = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc"))
    dist = dists[initial_time]
    ics = initial_conditions(dist, [1], 100, "sorted", ref_itp)

    tspan = (0.0, 100.0)

    cp_default = ClumpParameters(ref_itp) 
    cp = ClumpParameters(ref_itp, α, τ, cp_default.R, cp_default.f)

    sp = SpringParameters(k -> 0.0, 0.0)

    nw_type = "none"
    icons = form_connections(ics, nw_type)

    gdm = ImmortalModel()

    land = Land(verbose = false)

    rp = RaftParameters(
        ics = ics,
        clumps = cp,
        springs = sp,
        connections = icons,
        t0 = first(tspan),
        gd_model = gdm
    )

    prob = ODEProblem(WaterWind!, rp.ics, tspan, rp)

    sol = solve(
        prob, 
        Tsit5(), abstol = 1e-6, reltol = 1e-6,
        callback = CallbackSet(
            cb_update(showprogress = false), 
            callback(land), 
            callback(gdm), 
            cb_connections(network_type = nw_type))
        )

    return RaftTrajectory(sol, rp, ref_itp, dt = 0.1)
end

