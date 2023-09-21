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

# x_range, y_range = sph2xy(range(x0 - 1, x0 + 1, step = 0.5), range(y0 - 1, y0 + 1, step = 0.5), ref_itp)
x_range, y_range = sph2xy(range(x0 - 2, x0 + 4, step = 0.3), range(y0 - 4, y0 + 3, step = 0.3), ref_itp)
# x_range, y_range = sph2xy(range(x0 - 30, x0 + 5, step = 0.5), range(y0 - 7, y0 + 2, step = 0.5), ref_itp)
# tspan = (121.0, 122.0)
# tspan = (121.0, 151.0) # April 1 - May 1
tspan = (121.0, 212.0) # April 1 - July 1
cp = ClumpParameters(ref_itp)

# spring_k_constant = x -> 5

function spring_k(x::Real; A::Real = 5.0, k10::Real = 2*step(x_range))
    return A * (5/k10) * x * exp(1 - (5/k10)*x)
end

# gd_model = ImmortalModel()
gd_model = BrooksModel(params = BrooksModelParameters(temp_itp, no3_itp, clumps_limits = (0, 1000)), verbose = true)
# gd_model = BrooksModel(verbose = true)
rp = RaftParameters(x_range, y_range, cp, spring_k, first(tspan), "full", gd_model)

prob_raft = ODEProblem(Raft!, rp.ics, tspan, rp)

@info "Solving model."

land = Land(verbose = true)

@time sol_raft = solve(prob_raft, 
    Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(cb_loc2label(), callback(land), callback(gd_model))
);

rtr = RaftTrajectory(sol_raft, rp, ref_itp)

@info "Generating reference clump."

gd_model = ImmortalModel()
land = Land(verbose = true)
rp_1c = OneClumpRaftParameters(sph2xy(x0, y0, ref_itp)..., cp, first(tspan), gd_model)
prob_raft_1c = ODEProblem(Raft!, rp_1c.ics, tspan, rp_1c)

@time sol_clump = solve(prob_raft_1c, 
    Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(cb_loc2label(), callback(land), callback(gd_model))
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

