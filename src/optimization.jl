include(joinpath(@__DIR__, "SargassumBOMB.jl"))

using Surrogates
############################

isdefined(@__MODULE__, :OPTIMIZATION_PARAMETER_NAMES) || (const OPTIMIZATION_PARAMETER_NAMES = ["δ", "a", "β", "A_spring", "μ_max", "m", "k_N"])

"""
    mutable struct OptimizationParameter{T}
"""
mutable struct OptimizationParameter{T<:Real} 
    name::String
    default::T
    bounds::Tuple{T, T}
    val::T
    opt::Union{Nothing, T}
    optimizable::Bool

    function OptimizationParameter(name::String, default::Real, bounds::Tuple{Real, Real}, optimizable::Bool)
        @assert name in OPTIMIZATION_PARAMETER_NAMES
        @assert first(bounds) < last(bounds)

        def, lb, ub = promote(default, first(bounds), last(bounds))

        return new{typeof(def)}(name, def, (lb, ub), def, nothing, optimizable)
    end
end


"""
    mutable struct BOMBOptimizationProblem{T, U}
"""
mutable struct BOMBOptimizationProblem{T<:Real, U<:Integer}
    params::Dict{String, OptimizationParameter{T}}
    rhs::Function
    immortal::Bool
    tspan::Tuple{Tuple{U, U}, Tuple{U, U}}
    n_levels::U
    t_extra::U
    opt::Union{Nothing, T}
    seed::U

    function BOMBOptimizationProblem(;
        params::Vector{<:OptimizationParameter{T}},
        rhs::Function,
        immortal::Bool,
        tspan::Tuple{Tuple{U, U}, Tuple{U, U}},
        n_levels::U,
        t_extra::U,
        seed::U = 1234) where {T<:Real, U<:Integer}

        @assert length(params) > 0 "Must optimize at least one parameter."
        @assert rhs in [Raft!, WaterWind!]
        @assert DateTime(first(tspan)...) < DateTime(last(tspan)...)
        @assert allunique([param.name for param in params]) "Each parameter can only appear once."
        @assert n_levels > 0 "Need at least one level"

        params_dict = Dict(param.name => param for param in params)

        return new{T, U}(params_dict, rhs, immortal, tspan, n_levels, t_extra, nothing, seed)
    end
end

"""
    integrate_bomb(bop::BOMBOptimizationProblem; type::String)
"""
function integrate_bomb(bop::BOMBOptimizationProblem; type::String = "val")
    seed!(bop.seed)

    @assert type in ["val", "default", "opt"]
    if type == "opt"
        @assert bop.opt !== nothing "BOMBOptimizationProblem must be optimized to use optimal values."
    end

    initial_time = first(bop.tspan)
    final_time = last(bop.tspan)
    if type == "val"
        δ, a, β, A_spring, μ_max, m, k_N = [bop.params[param].val for param in OPTIMIZATION_PARAMETER_NAMES]
    elseif type == "default"
        δ, a, β, A_spring, μ_max, m, k_N = [bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES]
    elseif type == "opt"
        δ, a, β, A_spring, μ_max, m, k_N = [bop.params[param].optimizable ? bop.params[param].opt : bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES]
    end

    # TIME
    tstart = Day(DateTime(initial_time...) - DateTime(yearmonth(water_itp.time_start)...)).value |> float
    tend = tstart + Day(DateTime(final_time...) - DateTime(initial_time...)).value + bop.t_extra
    tspan = (tstart, tend)

    @assert tend > tstart "t_extra too negative"

    # ICS
    dist = DISTS_2018[initial_time]
    ics = initial_conditions(dist, [1], bop.n_levels, "levels", ref_itp)

    # CLUMPS
    cp = ClumpParameters(ref_itp, 
        δ = δ, 
        a = a,
        β = β)

    # SPRINGS
    p1 = sph2xy(dist.lon[1], dist.lat[1], ref_itp)
    p2 = sph2xy(dist.lon[2], dist.lat[2], ref_itp)
    ΔL = norm(p1 - p2)
    
    k10 = 2*ΔL
    L_spring = k10/5
    function spring_k(x::Real; A::Real = A_spring, k10::Real = k10)
        return A * (5/k10) * x * exp(1 - (5/k10)*x)
    end
    sp = SpringParameters(spring_k, L_spring)

    # CONNECTIONS
    nw_type = "nearest"
    n_conn = 10
    icons = form_connections(ics, nw_type, neighbor_parameter = n_conn)

    # BIOLOGY
    if bop.immortal
        gdm = ImmortalModel()
    else
        bmp = BrooksModelParameters(temp_itp, no3_itp, clumps_limits = (0, 10000), 
            μ_max = μ_max,
            m = m,
            k_N = k_N)
        gdm = BrooksModel(params = bmp)
    end

    land = Land()

    rp = RaftParameters(
        ics = ics,
        clumps = cp,
        springs = sp,
        connections = icons,
        t0 = first(tspan),
        gd_model = gdm
    )

    prob = ODEProblem(bop.rhs, rp.ics, tspan, rp)

    sol = solve(
        prob, 
        Tsit5(), abstol = 1e-6, reltol = 1e-6,
        callback = CallbackSet(
            cb_update(showprogress = false), 
            callback(land), 
            callback(gdm), 
            cb_connections(network_type = nw_type, neighbor_parameter = n_conn))
        )

    return (RaftTrajectory(sol, rp, ref_itp, dt = 0.1), tstart, tend)
end

"""
    loss_bomb(u, bop::BOMBOptimizationProblem)
"""
function loss_bomb(u, bop::BOMBOptimizationProblem)
    optimizable = [bop.params[param].name for param in keys(bop.params) if bop.params[param].optimizable]

    for i = 1:length(optimizable)
        bop.params[optimizable[i]].val = u[i]
    end

    initial_time = first(bop.tspan)
    final_time = last(bop.tspan)

    target = DISTS_2018[final_time].sargassum[:,:,1]
    target = target/sum(target)

    rtr, tstart, tend = integrate_bomb(bop)
    rtr = time_slice(rtr, (tend - bop.t_extra, tend))
    data = bins(rtr, DISTS_2018[final_time])
    data = data/sum(data)

    return sum((data - target) .^ 2)   
end

"""
    loss_bomb(bop::BOMBOptimizationProblem, type::String)
"""
function loss_bomb(bop::BOMBOptimizationProblem, type::String)

    initial_time = first(bop.tspan)
    final_time = last(bop.tspan)

    target = DISTS_2018[final_time].sargassum[:,:,1]
    target = target/sum(target)

    rtr, tstart, tend = integrate_bomb(bop, type = type)
    rtr = time_slice(rtr, (tend - bop.t_extra, tend))
    data = bins(rtr, DISTS_2018[final_time])
    data = data/sum(data)

    return sum((data - target) .^ 2)   
end

"""
    surrogate_bomb(n_samples_sur, bop::BOMBOptimizationProblem)
"""
function surrogate_bomb(n_samples_sur::Integer, bop::BOMBOptimizationProblem)
    lower_bound = typeof(bop).parameters[1][]
    upper_bound = typeof(bop).parameters[1][]

    optimizable = [bop.params[param].name for param in keys(bop.params) if bop.params[param].optimizable]

    for i = 1:length(optimizable)
        lb, ub = bop.params[optimizable[i]].bounds            
        push!(lower_bound, lb)
        push!(upper_bound, ub)
    end

    xys = Surrogates.sample(n_samples_sur, lower_bound, upper_bound, SobolSample())

    zs = zeros(eltype(xys).parameters[1], length(xys))

    Threads.@threads for i = 1:length(xys)
        zs[i] = loss_bomb(xys[i], bop)
    end

    rb = RadialBasis(xys, zs, lower_bound, upper_bound)

    return rb
end

"""
    optimize_bomb(u; bop::BOMBOptimizationProblem)
"""
function optimize_bomb(radial_basis::RadialBasis, bop::BOMBOptimizationProblem; maxiters_opt::Integer = 50)

    opt_vals = surrogate_optimize(u -> loss_bomb(u, bop), 
                        DYCORS(), 
                        radial_basis.lb, radial_basis.ub, radial_basis, SobolSample(), 
                        maxiters = maxiters_opt)

    optimizable = [bop.params[param].name for param in keys(bop.params) if bop.params[param].optimizable]

    for i = 1:length(optimizable)
        bop.params[optimizable[i]].opt = opt_vals[1][i]
    end

    bop.opt = opt_vals[2]

    return bop
end

"""
    plot_bop(bop::BOMBOptimizationProblem)
"""
function plot_bop(bop::BOMBOptimizationProblem)
    if bop.opt === nothing
        @warn "BOMBOptimizationProblem is not optimized; showing defaults."
    end

    initial_time = first(bop.tspan)
    final_time = last(bop.tspan)

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
    # initial distribution (SIMUL, unoptimized)
    ax = geo_axis(fig[2, 1], limits = limits, title = "SIMUL initial [default] $(monthname(initial_time[2])), week 1")
    rtr_dt, tstart, tend = integrate_bomb(bop, type = "default")
    dist = DISTS_2018[initial_time]
    rtr_dt_initial = time_slice(rtr_dt, (first(rtr_dt.t), first(rtr_dt.t)))
    trajectory_hist!(ax, rtr_dt_initial, dist)
    land!(ax)

    # final distribution (SIMUL, unoptimized)
    ax = geo_axis(fig[2, 2], limits = limits, title = "SIMUL final [default] $(monthname(final_time[2])), week 1")
    rtr_final = time_slice(rtr_dt, (tend - 8, tend))
    trajectory_hist!(ax, rtr_final, dist)
    land!(ax)

    ### OPTIMIZED
    # initial distribution (SIMUL, unoptimized)
    ax = geo_axis(fig[3, 1], limits = limits, title = "SIMUL initial [optim] $(monthname(initial_time[2])), week 1")
    rtr_dt, tstart, tend = integrate_bomb(bop, type = "opt")
    dist = DISTS_2018[initial_time]
    rtr_dt_initial = time_slice(rtr_dt, (first(rtr_dt.t), first(rtr_dt.t)))
    trajectory_hist!(ax, rtr_dt_initial, dist)
    land!(ax)

    # final distribution (SIMUL, unoptimized)
    ax = geo_axis(fig[3, 2], limits = limits, title = "SIMUL final [optim] $(monthname(final_time[2])), week 1")
    rtr_final = time_slice(rtr_dt, (tend - 8, tend))
    trajectory_hist!(ax, rtr_final, dist)
    land!(ax)

    # strings
    default_loss = loss_bomb(bop, "default")
    optimized_loss = loss_bomb(bop, "opt")

    ltx(x) = latexify(x, fmt = FancyNumberFormatter(4))

    dl_ltx, ol_ltx = ltx(default_loss), ltx(optimized_loss)
    ol_ltx = latexify(optimized_loss, fmt = FancyNumberFormatter(4))
    δ_def, a_def, β_def, A_spring_def, μ_max_def, m_def, k_N_def = ltx.([bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES])
    δ_opt, a_opt, β_opt, A_spring_opt, μ_max_opt, m_opt, k_N_opt = ltx.([bop.params[param].opt for param in OPTIMIZATION_PARAMETER_NAMES])

    if bop.rhs == WaterWind!
        fig[-3,:] = Label(fig, L"\text{WaterWind}")
    elseif bop.rhs == Raft!
        fig[-3,:] = Label(fig, L"\text{BOMB}")
    end

    fig[-2,:] = Label(fig, L"Loss(default) = %$(dl_ltx), Loss(opt) =  %$(ol_ltx)")

    fig[-1,:] = Label(fig, L"Defaults: $\delta =$ %$(δ_def), $a =$ %$(a_def), $\beta =$ %$(β_def), $A_\text{spring} =$ %$(A_spring_def), $\mu_\text{max} =$ %$(μ_max_def), $m =$ %$(m_def), $k_N =$ %$(k_N_def)")
    
    fig[0,:] = Label(fig, L"Optimals: $\delta =$ %$(δ_opt), $a =$ %$(a_opt), $\beta =$ %$(β_opt), $A_\text{spring} =$ %$(A_spring_opt), $\mu_\text{max} =$ %$(μ_max_opt), $m =$ %$(m_opt), $k_N =$ %$(k_N_opt)")

    outfile = joinpath(@__DIR__, "..", "figures", "opt_test.png")
    rm(outfile, force = true)
    save(outfile, fig)

    return fig
end

###############################################################################################

# ["δ", "a", "β", "A_spring", "μ_max", "m", "k_N"]

# OPTIMIZING
initial_time = (2018, 3)
final_time = (2018, 4)
t_extra = 7

δ_param = OptimizationParameter("δ",                1.25,   (1.05, 1.5),        false)
a_param = OptimizationParameter("a",                1.0e-4, (1.0e-5, 1.0e-3),   false)
β_param = OptimizationParameter("β",                0.0,    (0.0, 0.01),        true)
A_spring_param = OptimizationParameter("A_spring",  3.0,    (0.1, 10.0),        false)
μ_max_param = OptimizationParameter("μ_max",        0.1,    (0.05, 0.5),        false)
m_param = OptimizationParameter("m",                0.05,   (0.0, 0.1),         false)
k_N_param = OptimizationParameter("k_N",            0.012,  (0.005, 0.05),      false)

bop = BOMBOptimizationProblem(
    params = [δ_param, a_param, β_param, A_spring_param, μ_max_param, m_param, k_N_param],
    rhs = WaterWind!,
    immortal = true,
    tspan = (initial_time, final_time),
    n_levels = 2,
    t_extra = 7,
)

n_sur_samples = 50
maxiters_opt = 50

@info "Computing surrogate."
@time rb = surrogate_bomb(n_sur_samples, bop)
@info "Optimizing surrogate."
@time bop = optimize_bomb(rb, bop, maxiters_opt = maxiters_opt)
@info "Plotting."
@time plot_bop(bop)