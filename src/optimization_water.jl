include(joinpath(@__DIR__, "SargassumBOMB.jl"))

using Surrogates

########################################################################

function integrate_water(
    initial_time::NTuple{2, Integer},
    final_time::NTuple{2, Integer},
    t_extra::Real = 7;
    α::Real = ClumpParameters(ref_itp).α, 
    β::Real = ClumpParameters(ref_itp).β,
    seed::Integer = 1234)

    seed!(seed)

    tstart = Day(DateTime(initial_time...) - DateTime(yearmonth(water_itp.time_start)...)).value |> float
    tend = tstart + Day(DateTime(final_time...) - DateTime(initial_time...)).value + t_extra
    tspan = (tstart, tend)

    dist = DISTS_2018[initial_time]
    ics = initial_conditions(dist, [1], 2, "levels", ref_itp)

    cp_default = ClumpParameters(ref_itp) 
    cp = ClumpParameters(ref_itp, α, cp_default.τ, cp_default.R, cp_default.f, β)

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

    return (RaftTrajectory(sol, rp, ref_itp, dt = 0.1), tstart, tend)
end

function loss_water(
    initial_time::NTuple{2, Integer},
    final_time::NTuple{2, Integer},
    t_extra::Real = 7;
    α::Real = ClumpParameters(ref_itp).α, 
    β::Real = ClumpParameters(ref_itp).β,
    seed::Integer = 1234)

    target = DISTS_2018[final_time].sargassum[:,:,1]
    target = target/sum(target)

    rtr, tstart, tend = integrate_water(initial_time, final_time, t_extra, α = α, β = β, seed = seed)
    rtr = time_slice(rtr, (tend - 8, tend))
    data = bins(rtr, DISTS_2018[final_time])
    data = data/sum(data)

    return sum((data - target) .^ 2)
end

#################################################################
# OPTIMIZING
initial_time = (2018, 3)
final_time = (2018, 6)
t_extra = 7
loss_opt(u) = loss_water(initial_time, final_time, t_extra, α = u[1], β = u[2])

n_samples_sur = 100         # default 100
maxiters_opt = 50           # default 50
lower_bound = [0.0, 0.0]    # [α, β] lower
upper_bound = [0.05, 0.05]  # [α, β] upper

@info "Computing surrogate"
xys = Surrogates.sample(n_samples_sur, lower_bound, upper_bound, SobolSample())
@time zs = loss_opt.(xys)

radial_basis = RadialBasis(xys, zs, lower_bound, upper_bound)

@info "Optimizing surrogate"
@time sol_sur = surrogate_optimize(loss_opt, DYCORS(), lower_bound, upper_bound, 
                                    radial_basis, SobolSample(), maxiters = maxiters_opt)


α_opt, β_opt = sol_sur[1]
optimized_loss = sol_sur[2]
default_loss = loss_water(initial_time, final_time, t_extra)

@info "Optimal params: α = $(α_opt), β = $(β_opt)"
@info "Default loss: $(default_loss)"
@info "Optmzed loss: $(optimized_loss)"

#################################################################
### PLOTTING

fig = Figure(
    # resolution = (1920, 1080), 
    resolution = (2220, 2120),
    fontsize = 50,
    figure_padding = (5, 100, 5, 5))

limits = (-100, -40, 5, 35)

### AFAI
# initial distribution (AFAI)
ax = geo_axis(fig[1, 1], limits = limits, title = "AFAI initial $(monthname(initial_time[2])), week 1")
SFA_plot!(ax, initial_time, 1)
land!(ax)

# final distribution (AFAI)
ax = geo_axis(fig[1, 2], limits = limits, title = "AFAI final $(monthname(final_time[2])), week 1")
SFA_plot!(ax, final_time, 1)
land!(ax)

### UNOPTIMIZED
# initial distribution (WATER, unoptimized)
ax = geo_axis(fig[2, 1], limits = limits, title = "WATER initial [default] $(monthname(initial_time[2])), week 1")
rtr_dt, tstart, tend = integrate_water(initial_time, final_time, t_extra)
dist = DISTS_2018[initial_time]
rtr_dt_initial = time_slice(rtr_dt, (first(rtr_dt.t), first(rtr_dt.t)))
trajectory_hist!(ax, rtr_dt_initial, dist)
land!(ax)

# final distribution (WATER, unoptimized)
ax = geo_axis(fig[2, 2], limits = limits, title = "WATER final [default] $(monthname(final_time[2])), week 1")
rtr_final = time_slice(rtr_dt, (tend - 8, tend))
trajectory_hist!(ax, rtr_final, dist)
land!(ax)

### OPTIMIZED
# initial distribution (WATER, unoptimized)
ax = geo_axis(fig[3, 1], limits = limits, title = "WATER initial [optim] $(monthname(initial_time[2])), week 1")
rtr_dt, tstart, tend = integrate_water(initial_time, final_time, t_extra, 
                                        α = α_opt, 
                                        β = β_opt)
dist = DISTS_2018[initial_time]
rtr_dt_initial = time_slice(rtr_dt, (first(rtr_dt.t), first(rtr_dt.t)))
trajectory_hist!(ax, rtr_dt_initial, dist)
land!(ax)

# final distribution (WATER, unoptimized)
ax = geo_axis(fig[3, 2], limits = limits, title = "WATER final [optim] $(monthname(final_time[2])), week 1")
rtr_final = time_slice(rtr_dt, (tend - 8, tend))
trajectory_hist!(ax, rtr_final, dist)
land!(ax)

# strings
dl_ltx = latexify(default_loss, fmt = FancyNumberFormatter(4))
ol_ltx = latexify(optimized_loss, fmt = FancyNumberFormatter(4))
α_opt_ltx = latexify(α_opt, fmt = FancyNumberFormatter(4))
β_opt_ltx = latexify(β_opt, fmt = FancyNumberFormatter(4))

fig[-2,:] = Label(fig, L"\text{WATER}")
fig[-1,:] = Label(fig, L"$\alpha = $ %$(α_opt_ltx), $\beta =$ %$(β_opt_ltx)")
fig[0,:] = Label(fig, L"L Default = %$(dl_ltx), L Opt =  %$(ol_ltx)")

outfile = joinpath(@__DIR__, "..", "figures", "water_test.png")
rm(outfile, force = true)
save(outfile, fig)

fig