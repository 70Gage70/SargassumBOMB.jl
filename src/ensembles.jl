include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

using SargassumFromAFAI
using Random: seed!
######################################################################

function ensemble(
    start_date::NTuple{2, Int64}, 
    end_date::NTuple{2, Int64};
    rhs::Function = Raft!,
    cp::ClumpParameters = ClumpParameters(EQR), 
    cb_connections_type::String = "nearest", 
    rtr_dt::Real = 1.0)

    dists = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc"))
    dist = dists[start_date]

    tstart = Day(DateTime(start_date...) - DateTime(yearmonth(water_itp.time_start)...)).value |> float
    tend = tstart + Day(DateTime(end_date...) - DateTime(start_date...)).value
    tspan = (tstart, tend)

    @info "Integrating $(tspan)"

    # cp = ClumpParameters(EQR)
    # cp = ClumpParameters(EQR, δ = 2.0)

    ###################################################################### SPRINGS
    # x_range = range(-65.0, -55.0, step = 0.2)
    # y_range = range(8.0, 17.0, step = 0.2)
    # x_range, y_range = sph2xy(x_range, y_range, EQR)
    # ΔL = norm([x_range[1], y_range[1]] - [x_range[2], y_range[2]])

    p1 = sph2xy(dist.lon[1], dist.lat[1], EQR)
    p2 = sph2xy(dist.lon[2], dist.lat[2], EQR)
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
    # gdm = BrooksModel(params = BrooksModelParameters(TEMP_ITP.x, NO3_ITP.x, clumps_limits = (0, 1000)), verbose = true)
    # gdm = BrooksModel(verbose = true)

    ###################################################################### CONDITIONS

    # ics = initial_conditions(dist, [1], 100, "sorted", EQR)
    # ics = initial_conditions(dist, [1], 1, "uniform", EQR)
    ics = initial_conditions(dist, [1], 1000, "sample", EQR)

    # ics = initial_conditions(x_range, y_range)

    # icons = form_connections(ics, "nearest", neighbor_parameter = 10)
    icons = form_connections(ics, "radius", neighbor_parameter = k10)
    # icons = form_connections(ics, "full")
    # icons = form_connections(ics, "none")

    rp = RaftParameters(
        ics = ics,
        clumps = cp,
        springs = sp,
        connections = icons,
        t0 = first(tspan),
        gd_model = gdm
    )

    prob_raft = ODEProblem(rhs, rp.ics, tspan, rp)

    land = Land(verbose = false)

    @assert cb_connections_type in ["nearest", "none"]
    if cb_connections_type == "nearest"
        cb_c = cb_connections()
    elseif cb_connections_type == "none"
        cb_c = cb_connections(network_type = "none")
    end
        
    @time sol_raft = solve(
        prob_raft, 
        Tsit5(),
        # Tsit5(), abstol = 1e-6, reltol = 1e-6,
        # RK4(), dtmin = 0.1, dtmax = 0.1, force_dtmin = true,
        callback = CallbackSet(
            cb_update(showprogress = true), 
            callback(land), 
            callback(gdm), 
            cb_c)
        )

    return RaftTrajectory(sol_raft, rp, EQR, dt = rtr_dt)
end

dists = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc"))

april_plot = SargassumFromAFAI.plot(dists[(2018, 4)], size = (1920, 1080), legend = false, limits = (-100, -50, 5, 35))
july_plot = SargassumFromAFAI.plot(dists[(2018, 7)], size = (1920, 1080), legend = false, limits = (-100, -50, 5, 35))

# integrate to August 1, point is that the distribution for July take into account the entirety of 
# July, not just July 1st.
# rtrs = [ensemble((2018, t_start), (2018, 8)) for t_start = 4:6]
# rtrs = [ensemble(Leeway!, (2018, 4), (2018, 6))]
# rtrs = [ensemble(Raft!, (2018, 4), (2018, 6))]

# @info "Plotting results."

# limits = (-100, -50, 5, 35) # full plot``

# fig_tot = default_fig()
# ax = geo_axis(fig_tot[1, 1], limits = limits, title = L"\mathrm{April/May/June/July TOT - July}")

# ### hist

# lon_bins = range(-100, -50, length = 100)
# lat_bins = range(5, 35, length = 100)
# tr_hist = trajectory_hist!(ax, rtrs, lon_bins, lat_bins)

# land!(ax)

##################################################################################################

# t_start = last(rtrs[1].t) - 31.0
# t_end = last(rtrs[1].t)
# rtrs_slice = [time_slice(rtr, (t_start, t_end)) for rtr in rtrs]

# limits = (-100, -50, 5, 35) # full plot``

# fig_jul = default_fig()
# ax = geo_axis(fig_jul[1, 1], limits = limits, title = L"\mathrm{APRIL}")

# ### hist

# lon_bins = range(-100, -50, length = 100)
# lat_bins = range(5, 35, length = 100)
# tr_hist = trajectory_hist!(ax, rtrs_slice, lon_bins, lat_bins)

# land!(ax)

# fig_jul

#########

###

# cp_default = ClumpParameters(EQR)
cp_default = ClumpParameters(EQR, δ = 2.0)
cp_water = ClumpParameters(EQR, 0.0, 0.0, 0.0, 0.0, 0.0)
cp_wind = ClumpParameters(EQR, cp_default.α, 0.0, 0.0, 0.0, 0.0)

seed!(1234)
rtr_water = ensemble((2018, 4), (2018, 8), rhs = Leeway!, cp = cp_water, rtr_dt = 0.1)
seed!(1234)
rtr_wind = ensemble((2018, 4), (2018, 8), rhs = Leeway!, cp = cp_wind, rtr_dt = 0.1)
seed!(1234)
rtr_none = ensemble((2018, 4), (2018, 8), rhs = Raft!, cp = cp_default, cb_connections_type = "none", rtr_dt = 0.1)
seed!(1234)
rtr_near = ensemble((2018, 4), (2018, 8), rhs = Raft!, cp = cp_default, cb_connections_type = "nearest", rtr_dt = 0.1)

fig = default_fig()
limits = (-100, -50, 5, 35)


function change_trh(tspan)
    for i in 1:length(fig.content)
        delete!(fig.content[1])
    end

    day = round(first(tspan) - 90 + 1, digits = 2)
    lon_bins = range(-100, -50, length = 50)
    lat_bins = range(5, 35, length = 50)

    ax_water = geo_axis(fig[1, 1], limits = limits, title = L"\mathrm{Water: Day } \, %$(day)")
    trajectory_hist!(ax_water, time_slice(rtr_water, tspan), lon_bins, lat_bins, opts= (
        colormap = Reverse(:RdYlGn),
        colorscale = x -> x == 0.0 ? -1.0 : x))
    land!(ax_water)

    ax_wind = geo_axis(fig[1, 2], limits = limits, title = L"\mathrm{Water + Wind: Day } \, %$(day)")
    trajectory_hist!(ax_wind, time_slice(rtr_wind, tspan), lon_bins, lat_bins, opts= (
        colormap = Reverse(:RdYlGn),
        colorscale = x -> x == 0.0 ? -1.0 : x))
    land!(ax_wind)    

    ax_none = geo_axis(fig[2, 1], limits = limits, title = L"\mathrm{BOM: Day } \, %$(day)")
    trajectory_hist!(ax_none, time_slice(rtr_none, tspan), lon_bins, lat_bins, opts= (
        colormap = Reverse(:RdYlGn),
        colorscale = x -> x == 0.0 ? -1.0 : x))
    land!(ax_none) 
    
    ax_near = geo_axis(fig[2, 2], limits = limits, title = L"\mathrm{eBOM: Day } \, %$(day)")
    trajectory_hist!(ax_near, time_slice(rtr_near, tspan), lon_bins, lat_bins, opts= (
        colormap = Reverse(:RdYlGn),
        colorscale = x -> x == 0.0 ? -1.0 : x))
    land!(ax_near)     
end

trh_iterator = [(90 + 0.1*i, 90 + 0.1*(i + 1)) for i = 0:10:1219]
# trh_iterator = [(90 + 0.1*i, 90 + 0.1*(i + 1)) for i = 0:0]

# @time record(change_trh, fig, joinpath(@__DIR__, "..", "figures", "comparison-test-delta2.mp4"), trh_iterator; framerate = 15)

###

fig = default_fig()
limits = (-100, -50, 5, 35)

ax_water = geo_axis(fig[1, 1], limits = limits, title = L"\mathrm{Water}")
trajectory!(ax_water, rtr_water)
land!(ax_water)

ax_wind = geo_axis(fig[1, 2], limits = limits, title = L"\mathrm{Water + Wind}")
trajectory!(ax_wind, rtr_wind)
land!(ax_wind)    

ax_none = geo_axis(fig[2, 1], limits = limits, title = L"\mathrm{BOM}")
trajectory!(ax_none, rtr_none)
land!(ax_none) 

ax_near = geo_axis(fig[2, 2], limits = limits, title = L"\mathrm{eBOM}")
trajectory!(ax_near, rtr_near)
land!(ax_near)

supertitle = Label(fig[0, :], L"\delta = 2.0", fontsize = 50)

fig