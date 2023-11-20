function simulate(rp::RaftParameters)
    prob_raft = ODEProblem(Raft!, rp.ics.ics, tspan, rp)

    land = Land(verbose = false)

    sol_raft = solve(
        prob_raft, 
        Tsit5(),
        # Tsit5(), abstol = 1e-6, reltol = 1e-6,
        callback = CallbackSet(
            cb_update(showprogress = true), 
            callback(land), 
            callback(gdm), 
            cb_connections())
        )

@info "Generating model."

start_date = (2018, 4)
# end_date = (2018, 7)
end_date = (2018, 8)

DISTS_2018 = SargassumDistribution(SargassumFromAFAI.EXAMPLE_DIST_2018)
dist = DISTS_2018[start_date]

tstart = Day(DateTime(start_date...) - DateTime(yearmonth(WATER_ITP.x.time_start)...)).value |> float
tend = tstart + Day(DateTime(end_date...) - DateTime(start_date...)).value
tspan = (tstart, tend)

@info "Integrating $(tspan)"

cp = ClumpParameters(EQR_DEFAULT)
# cp = ClumpParameters(EQR_DEFAULT, δ = 3.0)

###################################################################### SPRINGS
# x_range = range(-65.0, -55.0, step = 0.2)
# y_range = range(8.0, 17.0, step = 0.2)
# y_range = range(15.0, 17.0, step = 0.2)
# x_range, y_range = sph2xy(x_range, y_range, EQR_DEFAULT)
# ΔL = norm([x_range[1], y_range[1]] - [x_range[2], y_range[2]])

p1 = sph2xy(dist.lon[1], dist.lat[1], EQR_DEFAULT)
p2 = sph2xy(dist.lon[2], dist.lat[2], EQR_DEFAULT)
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
# gdm = BrooksModel(params = BrooksModelParameters(TEMP_ITP, NO3_ITP, clumps_limits = (0, 1000)), verbose = true)
# gdm = BrooksModel(verbose = true)

###################################################################### CONDITIONS
seed!(1234)
# ics = InitialConditions(dist, [1], 1000, "sorted", EQR_DEFAULT)
# ics = InitialConditions(dist, [1], 1, "uniform", EQR_DEFAULT)
# ics = InitialConditions(dist, [1], 1000, "sample", EQR_DEFAULT)
ics = InitialConditions(dist, [1], 10, "levels", EQR_DEFAULT)

# ics = InitialConditions(x_range, y_range)

# icons = ConnectionsNearest(10)
# icons = ConnectionsRadius(k10)
icons = ConnectionsFull()
# icons = ConnectionsNone()

rp = RaftParameters(
    ics = ics,
    clumps = cp,
    springs = sp,
    connections = icons,
    t0 = first(tspan),
    gd_model = gdm
)

prob_raft = ODEProblem(Raft!, rp.ics.ics, tspan, rp)

land = Land(verbose = false)
    
@time sol_raft = solve(
    prob_raft, 
    Tsit5(),
    # Tsit5(), abstol = 1e-6, reltol = 1e-6,
    callback = CallbackSet(
        cb_update(showprogress = true), 
        callback(land), 
        cb_growth_death(), 
        cb_connections())
    )

rtr = RaftTrajectory(sol_raft, rp, EQR_DEFAULT, dt = 0.1)


