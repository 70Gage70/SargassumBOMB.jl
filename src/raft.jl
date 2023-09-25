include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

######################################################################

@info "Generating model."

# x0, y0 = -90, 23 # GoM
# x0, y0 = -64, 14
x0, y0 = -60, 13

Δ_clump = 0.3
x_range, y_range = sph2xy(range(x0 - 2, x0 + 4, step = Δ_clump), range(y0 - 4, y0 + 3, step = Δ_clump), ref_itp)

# tspan = (121.0, 122.0)
# tspan = (121.0, 151.0) # April 1 - May 1
tspan = (121.0, 212.0) # April 1 - July 1

cp = ClumpParameters(ref_itp)

# spring_k_constant = x -> 5
# L_spring = 2*(step(x_range) + step(y_range))/2
# sp = SpringParameters(spring_k_constant, L_spring)

k10 = 2*(step(x_range) + step(y_range))/2
L_spring = k10/5
function spring_k(x::Real; A::Real = 3.0, k10::Real = k10)
    return A * (5/k10) * x * exp(1 - (5/k10)*x)
end
sp = SpringParameters(spring_k, L_spring)

gdm = ImmortalModel()
# gdm = BrooksModel(params = BrooksModelParameters(temp_itp, no3_itp, clumps_limits = (0, 1000)), verbose = true)
# gdm = BrooksModel(verbose = true)

# rp = RectangularRaftParameters(x_range, y_range, cp, spring_k, first(tspan), "full", gdm)

ics = initial_conditions(x_range, y_range)
# icons = initial_connections(ics, "nearest", neighbor_parameter = 4)
# icons = initial_connections(ics, "full")
icons = initial_connections(ics, "none")
rp = RaftParameters(
    ics = ics,
    clumps = cp,
    springs = sp,
    connections = icons,
    t0 = first(tspan),
    gd_model = gdm
)

prob_raft = ODEProblem(Raft!, rp.ics, tspan, rp)

# wp = RaftParameters(x_range, y_range, cp, spring_k, first(tspan), "full", ImmortalModel())
# prob_raft = ODEProblem(Water!, rp.ics, tspan, rp)

@info "Solving model."

land = Land(verbose = true)

@time sol_raft = solve(prob_raft, 
    Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(cb_loc2label(), callback(land), callback(gdm))
);

rtr = RaftTrajectory(sol_raft, rp, ref_itp)

@info "Generating reference clump."

gdm = ImmortalModel()
land = Land(verbose = true)

ics_1c = initial_conditions(x0, y0, ref = ref_itp)
icons_1c = initial_connections(ics, "none")
rp_1c = RaftParameters(
    ics = ics_1c,
    clumps = cp,
    springs = sp,
    connections = icons_1c,
    t0 = first(tspan),
    gd_model = gdm
)

prob_raft_1c = ODEProblem(Raft!, rp_1c.ics, tspan, rp_1c)

@time sol_clump = solve(prob_raft_1c, 
    Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(cb_loc2label(), callback(land), callback(gdm))
);

ctr = RaftTrajectory(sol_clump, rp_1c, ref_itp)

@info "Plotting results."

limits = (-100, -50, 5, 35)
# limits = (-70, -50, 5, 25)

fig_COM = default_fig()
ax = geo_axis(fig_COM[1, 1], limits = limits, title = L"\mathrm{Raft COM}")

### raft
# trajectory!(ax, rtr)

### COM
# trajectory!(ax, rtr.com, 
#     opts = (linestyle = :dot, color = rtr.com.t, colormap = :heat, linewidth = 5)
# ) 

### COM 5
# rtr_5 = RaftTrajectory(sol_raft, rp, ref_itp, dt = 5.0)

# trajectory!(ax, rtr_5.com, 
#     opts = (linestyle = :dot, color = rtr_5.com.t, colormap = :viridis, linewidth = 5)
# ) 

# scatter!(ax, rtr_5.com.xy[:,1], rtr_5.com.xy[:,2])

### clump
# trajectory!(ax, ctr, 
#     opts = (color = ctr.t, colormap = :heat, linewidth = 2)
# ) 

### hist
rtr_dt = RaftTrajectory(sol_raft, rp, ref_itp, dt = 1.0)

lon_bins = range(-100, -50, length = 100)
lat_bins = range(5, 35, length = 100)
tr_hist = trajectory_hist!(ax, rtr_dt, lon_bins, lat_bins)


land!(ax)

### days legend
# tticks = collect(range(start = minimum(rtr.t), stop = maximum(rtr.t), length = 5))
# data_legend!(fig_COM[1,2], L"\mathrm{Days}", ticks = tticks)

### counts legend
min_cts, max_cts = getindex(tr_hist.colorrange)
tticks = collect(range(start = 0.0, stop = log10(max_cts), length = 5))
data_legend!(fig_COM[1,2], L"\log_{10} \left(\mathrm{Counts}\right)", ticks = tticks, colormap = Reverse(:RdYlGn))

fig_COM

