include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

using SargassumFromAFAI
######################################################################

@info "Generating model."

dists = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc"))
dist = dists[(2018, 4)]

# tspan = (121.0, 121.1)
tspan = (121.0, 212.0) # April 1 - July 1

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

# icons = form_connections(ics, "nearest", neighbor_parameter = 4)
icons = form_connections(ics, "full")
# icons = form_connections(ics, "none")

rp = RaftParameters(
    ics = ics,
    clumps = cp,
    springs = sp,
    connections = icons,
    t0 = first(tspan),
    gd_model = gdm
)

prob_raft = ODEProblem(Raft!, rp.ics, tspan, rp)

@info "Solving model."

land = Land(verbose = true)

@time sol_raft = solve(
    prob_raft, 
    Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(
        cb_update(showprogress = true), 
        callback(land), 
        callback(gdm), 
        cb_connections(radius = 2*k10))
);

rtr = RaftTrajectory(sol_raft, rp, ref_itp)

@info "Plotting results."

limits = (-100, -50, 5, 35) # full plot``
# limits = (-90, -38, 5, 22) # compare with dist

fig_COM = default_fig()
ax = geo_axis(fig_COM[1, 1], limits = limits, title = L"\mathrm{Raft COM}")

### hist
rtr_dt = RaftTrajectory(sol_raft, rp, ref_itp, dt = 1.0)

lon_bins = range(-100, -50, length = 100)
lat_bins = range(5, 35, length = 100)
tr_hist = trajectory_hist!(ax, rtr_dt, lon_bins, lat_bins)

land!(ax)

### counts legend
# min_cts, max_cts = getindex(tr_hist.colorrange)
# tticks = collect(range(start = 0.0, stop = log10(max_cts), length = 5))
# data_legend!(fig_COM[1,2], L"\log_{10} \left(\mathrm{Counts}\right)", ticks = tticks, colormap = Reverse(:RdYlGn))

fig_COM

