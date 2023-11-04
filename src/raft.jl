include(joinpath(@__DIR__, "SargassumBOMB.jl"))

######################################################################

@info "Generating model."

start_date = (2018, 4)
# end_date = (2018, 7)
end_date = (2018, 8)

dist = DISTS_2018[start_date]

tstart = Day(DateTime(start_date...) - DateTime(yearmonth(water_itp.time_start)...)).value |> float
tend = tstart + Day(DateTime(end_date...) - DateTime(start_date...)).value
tspan = (tstart, tend)

@info "Integrating $(tspan)"

cp = ClumpParameters(ref_itp)
# cp = ClumpParameters(ref_itp, δ = 3.0)

###################################################################### SPRINGS
# x_range = range(-65.0, -55.0, step = 0.2)
# y_range = range(8.0, 17.0, step = 0.2)
# y_range = range(15.0, 17.0, step = 0.2)
# x_range, y_range = sph2xy(x_range, y_range, ref_itp)
# ΔL = norm([x_range[1], y_range[1]] - [x_range[2], y_range[2]])

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
seed!(1234)
# ics = initial_conditions(dist, [1], 1000, "sorted", ref_itp)
# ics = initial_conditions(dist, [1], 1, "uniform", ref_itp)
# ics = initial_conditions(dist, [1], 1000, "sample", ref_itp)
ics = initial_conditions(dist, [1], 10, "levels", ref_itp)

# ics = initial_conditions(x_range, y_range)

# icons = form_connections(ics, "nearest", neighbor_parameter = 10)
# icons = form_connections(ics, "radius", neighbor_parameter = k10)
# icons = form_connections(ics, "full")
icons = form_connections(ics, "none")

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

# cb_c = cb_connections()
cb_c = cb_connections(network_type = "none")
    
@time sol_raft = solve(
    prob_raft, 
    Tsit5(),
    # Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(
        cb_update(showprogress = true), 
        callback(land), 
        callback(gdm), 
        cb_c)
    )

rtr = RaftTrajectory(sol_raft, rp, ref_itp, dt = 0.1)

# @info "Generating reference clump."

# gdm = ImmortalModel()
# land = Land(verbose = true)

# ics_1c = initial_conditions(x0, y0, ref = ref_itp)
# icons_1c = form_connections(ics, "none")
# rp_1c = RaftParameters(
#     ics = ics_1c,
#     clumps = cp,
#     springs = sp,
#     connections = icons_1c,
#     t0 = first(tspan),
#     gd_model = gdm
# )

# prob_raft_1c = ODEProblem(Raft!, rp_1c.ics, tspan, rp_1c)

# @time sol_clump = solve(prob_raft_1c, 
#     Tsit5(), abstol = 1e-6, reltol = 1e-6,
#     callback = CallbackSet(cb_loc2label(), callback(land), callback(gdm))
# );

# ctr = RaftTrajectory(sol_clump, rp_1c, ref_itp)

@info "Plotting results."

limits = (-100, -40, 5, 35)
# limits = (-70, -50, 5, 25)

fig_COM = default_fig()
ax = geo_axis(fig_COM[1, 1], limits = limits, title = L"\mathrm{Raft: April Week 1}")

### raft
# trajectory!(ax, rtr)

### COM
# trajectory!(ax, rtr.com, 
#     opts = (linestyle = :dot, color = rtr.com.t, colormap = :heat, linewidth = 5)
# ) 

### clump
# trajectory!(ax, ctr, 
#     opts = (color = ctr.t, colormap = :heat, linewidth = 2)
# ) 

### hist
rtr_dt = RaftTrajectory(sol_raft, rp, ref_itp, dt = 1.0)

dist = DISTS_2018[(2018, 4)]
rtr_dt_initial = time_slice(rtr_dt, (first(rtr_dt.t), first(rtr_dt.t)))
trajectory_hist!(ax, rtr_dt_initial, dist)

land!(ax)

### days legend
# tticks = collect(range(start = minimum(rtr.t), stop = maximum(rtr.t), length = 5))
# data_legend!(fig_COM[1,2], L"\mathrm{Days}", ticks = tticks)

### counts legend
# min_cts, max_cts = getindex(tr_hist.colorrange)
# tticks = collect(range(start = 0.0, stop = log10(max_cts), length = 5))
# data_legend!(fig_COM[1,2], L"\log_{10} \left(\mathrm{Counts}\right)", ticks = tticks, colormap = Reverse(:RdYlGn))

fig_COM

