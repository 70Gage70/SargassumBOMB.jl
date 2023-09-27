include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

using SargassumFromAFAI
######################################################################

function ensemble(start_date::NTuple{2, Int64}, end_date::NTuple{2, Int64}; rtr_dt::Real = 1.0)

    dists = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc"))
    dist = dists[start_date]

    tstart = Day(DateTime(start_date...) - DateTime(yearmonth(water_itp.time_start)...)).value |> float
    tend = tstart + Day(DateTime(end_date...) - DateTime(start_date...)).value
    tspan = (tstart, tend)

    @info "Integrating $(tspan)"

    cp = ClumpParameters(ref_itp)

    ###################################################################### SPRINGS
    p1 = sph2xy(dist.lon[1], dist.lat[1], ref_itp)
    p2 = sph2xy(dist.lon[2], dist.lat[2], ref_itp)
    ΔL = norm(p1 - p2)

    # spring_k_constant = x -> 5
    # sp = SpringParameters(spring_k_constant, ΔL)

    A_spring = 3.0
    k10 = 2*ΔL
    L_spring = k10/5
    function spring_k(x::Real; A::Real = A_spring, k10::Real = k10)
        return A * (5/k10) * x * exp(1 - (5/k10)*x)
    end
    sp = SpringParameters(spring_k, L_spring)

    ###################################################################### BIOLOGY

    gdm = ImmortalModel()
    # gdm = BrooksModel(params = BrooksModelParameters(temp_itp, no3_itp, clumps_limits = (0, 1000)), verbose = true)
    # gdm = BrooksModel(verbose = true)

    ###################################################################### CONDITIONS

    # ics = initial_conditions(dist, 100, "sorted", ref_itp)
    ics = initial_conditions(dist, 1, "uniform", ref_itp)

    # icons = initial_connections(ics, "nearest", neighbor_parameter = 4)
    icons = initial_connections(ics, "full")
    # icons = initial_connections(ics, "none")

    rp = RaftParameters(
        ics = ics,
        clumps = cp,
        springs = sp,
        connections = icons,
        t0 = first(tspan),
        gd_model = gdm
    )

    prob_raft = ODEProblem(Raft!, rp.ics, tspan, rp)

    land = Land(verbose = false)

    @time sol_raft = solve(
        prob_raft, 
        Tsit5(), abstol = 1e-6, reltol = 1e-6,
        callback = CallbackSet(
            cb_update(showprogress = true), 
            callback(land), 
            callback(gdm), 
            cb_connections_radius(radius = 2*k10))
    );

    return RaftTrajectory(sol_raft, rp, ref_itp, dt = rtr_dt)
end

# integrate to August 1, point is that the distribution for July take into account the entirety of 
# July, not just July 1st.
rtrs = [ensemble((2018, t_start), (2018, 8)) for t_start = 4:6]

@info "Plotting results."

limits = (-100, -50, 5, 35) # full plot``

fig_COM = default_fig()
ax = geo_axis(fig_COM[1, 1], limits = limits, title = L"\mathrm{April/May/June/July TOT - July}")

### hist

lon_bins = range(-100, -50, length = 100)
lat_bins = range(5, 35, length = 100)
tr_hist = trajectory_hist!(ax, rtrs, lon_bins, lat_bins)

land!(ax)

fig_COM

##################################################################################################

t_start = last(rtrs[1].t) - 31.0
t_end = last(rtrs[1].t)
rtrs_slice = [time_slice(rtr, (t_start, t_end)) for rtr in rtrs]

limits = (-100, -50, 5, 35) # full plot``

fig_COM = default_fig()
ax = geo_axis(fig_COM[1, 1], limits = limits, title = L"\mathrm{April/May/June/July JUL - July}")

### hist

lon_bins = range(-100, -50, length = 100)
lat_bins = range(5, 35, length = 100)
tr_hist = trajectory_hist!(ax, rtrs_slice, lon_bins, lat_bins)

land!(ax)

fig_COM