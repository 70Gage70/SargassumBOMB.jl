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
    integrate_bomb(bop::BOMBOptimizationProblem)
"""
function integrate_bomb(bop::BOMBOptimizationProblem)
    seed!(bop.seed)

    initial_time = first(bop.tspan)
    final_time = last(bop.tspan)

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
        δ = bop.params["δ"].val, 
        a = bop.params["a"].val,
        β = bop.params["β"].val)

    # SPRINGS
    p1 = sph2xy(dist.lon[1], dist.lat[1], ref_itp)
    p2 = sph2xy(dist.lon[2], dist.lat[2], ref_itp)
    ΔL = norm(p1 - p2)
    
    k10 = 2*ΔL
    L_spring = k10/5
    function spring_k(x::Real; A::Real = bop.params["A_spring"].val, k10::Real = k10)
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
            μ_max = bop.params["μ_max"].val,
            m = bop.params["m"].val,
            k_N = bop.params["k_N"].val)
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
    zs = [loss_bomb(xysi, bop) for xysi in xys]
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

###############################################################################################

# ["δ", "a", "β", "A_spring", "μ_max", "m", "k_N"]

# OPTIMIZING
initial_time = (2018, 3)
final_time = (2018, 6)
t_extra = 7

δ_param = OptimizationParameter("δ", 1.25, (1.05, 1.5), true)
a_param = OptimizationParameter("a", 1.0e-4, (1.0e-5, 1.0e-3), true)
β_param = OptimizationParameter("β", 0.0, (0.0, 0.01), true)
A_spring_param = OptimizationParameter("A_spring", 3.0, (0.1, 10.0), true)
μ_max_param = OptimizationParameter("μ_max", 0.1, (0.05, 0.5), true)
m_param = OptimizationParameter("m", 0.05, (0.0, 0.1), true)
k_N_param = OptimizationParameter("k_N", 0.012, (0.005, 0.05), true)

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
@time rb = surrogate_bomb(n_sur_samples, bop)
@time bop = optimize_bomb(rb, bop, maxiters_opt = maxiters_opt)