"""
    const OPTIMIZATION_PARAMETER_NAMES

A `Vector` of `String`s giving the names of all the parameters it is possible to optimize by default.

Equal to `["δ", "a", "σ", "A_spring", "μ_max", "m", "k_N"]`.
"""
const OPTIMIZATION_PARAMETER_NAMES = ["δ", "a", "σ", "A_spring", "μ_max", "m", "k_N"]

"""
    struct LossFunction

A container for a function used for measuring the distance between two binned histograms for \
optimization. 

### Fields 

- `f`: A `Function`. This function must be callable as `f(a, b)`, where `a` and `b` are matrices and it \
must return a `Real` such that `f(a, a) = 0.0`. Example: `(a, b) -> sum((a - b) .^ 2)`.
- `name`: A `String` giving the name of the loss function, e.g. `"L1"`.
"""
struct LossFunction
    f::Function
    name::String

    function LossFunction(; f::Function, name::String)

        m1_test = rand(4, 4)
        m2_test = rand(4, 4)

        try 
            f(m1_test, m2_test)
        catch e
            @warn "Could not evaluate `f(::Matrix, ::Matrix)`. Ensure that `f` can accept two
                    matrix arguments."
            throw(e)
        end
    
        @assert f(m1_test, m2_test) isa Real "The loss function must evaluate to a real number."
        @assert f(m1_test, m1_test) <= f(m1_test, m2_test) "`f` doesn't decrease when data are closer together."

        return new(f, name)
    end
end

"""
    const LOSS_L1

A `LossFunction` for the L1 norm.
"""
const LOSS_L1 = LossFunction(f = (a::Matrix, b::Matrix) -> sum((a - b) .^ 2) , name = "L1")

"""
    const LOSS_COR

A `LossFunction` for the negative correlation.
"""
const LOSS_COR = LossFunction(f = (a::Matrix, b::Matrix) -> -cor(vec(a), vec(b)), name = "-COR")

"""
    mutable struct OptimizationParameter{T}

A container for the data, values and bounds for a parameter to be optimized.

### Fields 

- `name`: A `String` with the name of the parameter, must be one of [`OPTIMIZATION_PARAMETER_NAMES`](@ref).
- `default`: The default (starting) value of the parameter.
- `bounds`: A `Tuple` giving the upper and lower bounds of the parameter.
- `val`: A current value of the parameter.
- `opt`: A optimal value of the parameter.
- `optimizable`: A `Bool` such that, if `false`, the parameter will not be optimized away from its default value.

### Constructor

Use 

`OptimizationParameter(name, default, bounds, optimizable; val = nothing, opt = nothing)`
"""
mutable struct OptimizationParameter{T<:Real} 
    name::String
    default::T
    bounds::Tuple{T, T}
    val::T
    opt::Union{Nothing, T}
    optimizable::Bool

    function OptimizationParameter(
        name::String, 
        default::Real, 
        bounds::Tuple{Real, Real},
        optimizable::Bool;
        val::Union{Nothing, Real} = nothing,
        opt::Union{Nothing, Real} = nothing)

        @assert name in OPTIMIZATION_PARAMETER_NAMES
        @assert first(bounds) < last(bounds)

        def, lb, ub = promote(default, first(bounds), last(bounds))

        op_val = val === nothing ? def : val
        op_opt = opt === nothing ? nothing : promote(opt, def)[1]

        return new{typeof(def)}(name, def, (lb, ub), op_val, op_opt, optimizable)
    end
end

"""
    mutable struct BOMBOptimizationProblem{T, U}

A container for all the data defining an optimization problem.

### Fields 

- `params`: A `Dict` mapping each element of `[OPTIMIZATION_PARAMETER_NAMES](@ref)` to an [`OptimizationParameter`](@ref) \
that contains it.
- `rhs`: The `Function` to integrate, generally should be [`Raft!`](@ref), but [`WaterWind!`] can be used \
for testing purposes.
- `immortal`: A `Bool` such that if `true`, the [`ImmortalModel`](@ref) will be used, resulting in no clump \
growths or deaths.
- `tspan`: A `Tuple` of the form `((year1, month1), (year2, month2))` giving the integration time span.
- `n_levels`: The number of levels to be used in the [`InitialConditions`](@ref) for the `SargassumDistribution`. \
The more levels, the more clumps initially exist.
- `t_extra`: A number of extra days to add to the integration at the end.
- `loss_func`: The [`LossFunction`](@ref) used during the optimization.
- `opt`: The minimal [`LossFunction`](@ref) obtained. If `nothing`, the problem is considered unoptimized.
- `seed`: A `Random.seed!` used in the integration.

### Constructor

Use 

`BOMBOptimizationProblem(; kwargs...)` where each field has a named kwarg.

If not provided, `loss_func` defaults to [`LOSS_COV`](@ref) and `seed` defaults to `1234`. In general, 
`opt` should not be provided directly but it can be to bypass certain checks.
"""
mutable struct BOMBOptimizationProblem{T<:Real, U<:Integer}
    params::Dict{String, OptimizationParameter{T}}
    rhs::Function
    immortal::Bool
    tspan::Tuple{Tuple{U, U}, Tuple{U, U}}
    n_levels::U
    t_extra::U
    loss_func::LossFunction
    opt::Union{Nothing, T}
    seed::U

    function BOMBOptimizationProblem(;
        params::Dict{String, OptimizationParameter{T}},
        rhs::Function,
        immortal::Bool,
        tspan::Tuple{Tuple{U, U}, Tuple{U, U}},
        n_levels::U,
        t_extra::U,
        loss_func::LossFunction = LOSS_COV,
        opt::Union{Nothing, T} = nothing,
        seed::U = 1234) where {T<:Real, U<:Integer}

        @assert length(params) > 0 "Must optimize at least one parameter."
        @assert rhs in [Raft!, WaterWind!]
        @assert DateTime(first(tspan)...) < DateTime(last(tspan)...)
        @assert n_levels > 0 "Need at least one level"

        return new{T, U}(params, rhs, immortal, tspan, n_levels, t_extra, loss_func, opt, seed)
    end
end

"""
    simulate(bop::BOMBOptimizationProblem)

Integrate `bop` by constructing the [`RaftParameters`](@ref) implied by its fields using [`simulate(::RaftParameters)`](@ref).
"""
function simulate(bop::BOMBOptimizationProblem)
    δ, a, σ, A_spring, μ_max, m, k_N = [bop.params[param].val for param in OPTIMIZATION_PARAMETER_NAMES]

    # time
    start_date, end_date = bop.tspan
    dist = SargassumDistribution(SargassumFromAFAI.EXAMPLE_DIST_2018)[start_date]
    tspan = yearmonth2tspan(start_date, end_date, t_extra = (0, bop.t_extra))
    
    # clumps
    clumps = ClumpParameters(δ = δ, a = a, σ = σ)
    
    # springs
    L_spring = ΔL(dist)
    function spring_k(x::Real; A::Real = A_spring, L::Real = L_spring)
        return A_spring * (exp((x - 2*L)/0.2) + 1)^(-1)
    end

    springs = SpringParameters(spring_k, L_spring)
    
    # growth-death
    if bop.immortal
        gd_model = ImmortalModel()
    else
        bmp = BrooksModelParameters(TEMP_ITP.x, NO3_ITP.x, 
            μ_max = μ_max,
            m = m,
            k_N = k_N)
        gd_model = BrooksModel(params = bmp)
    end
    
    # initial conditions
    ics = InitialConditions(dist, [1], bop.n_levels, "levels", EQR_DEFAULT)
    
    # connections
    connections = ConnectionsNearest(10)
    
    # land
    land = Land()
    
    # integrating
    
    rp = RaftParameters(
        tspan = tspan,
        ics = ics,
        clumps = clumps,
        springs = springs,
        connections = connections,
        gd_model = gd_model,
        land = land
    )
    
    return simulate(rp, rhs = bop.rhs, showprogress = false)
end


"""
    loss(u::Vector{<:Real}, bop::BOMBOptimizationProblem)

Compute the loss associated with parameter values `u` in `bop`.

The vector `u` should be the same length as the number of `bop.params` which are optimizable. Then, 
the `bop` is integrated with the values of each parameter set equal to the entries of `u` in the 
order defined by [`OPTIMIZATION_PARAMETER_NAMES`](@ref).
"""
function loss(u::Vector{<:Real}, bop::BOMBOptimizationProblem)
    optimizable = [bop.params[param].name for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]

    @assert length(u) == length(optimizable) "The vector u must have the same number of entries as the number of optimizable parameters."

    for i = 1:length(optimizable)
        bop.params[optimizable[i]].val = u[i]
    end

    start_date, end_date = bop.tspan
    tstart, tend = yearmonth2tspan(start_date, end_date, t_extra = (0, bop.t_extra))

    target_dist = SargassumDistribution(SargassumFromAFAI.EXAMPLE_DIST_2018)[end_date]
    target = target_dist.sargassum[:,:,1]
    target = target/sum(target)

    rtr = simulate(bop)
    rtr = time_slice(rtr, (tend - bop.t_extra, tend))
    data = bins(rtr, target_dist)
    data = data/sum(data)

    return bop.loss_func.f(data, target)
end

"""
    surrogate_bomb(n_samples_sur, bop::BOMBOptimizationProblem)
"""
function surrogate_bomb(n_samples_sur::Integer, bop::BOMBOptimizationProblem)
    lower_bound = typeof(bop).parameters[1][]
    upper_bound = typeof(bop).parameters[1][]

    optimizable = [bop.params[param].name for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]

    for i = 1:length(optimizable)
        lb, ub = bop.params[optimizable[i]].bounds            
        push!(lower_bound, lb)
        push!(upper_bound, ub)
    end

    xys = Surrogates.sample(n_samples_sur, lower_bound, upper_bound, SobolSample())

    zs = zeros(eltype(xys).parameters[1], length(xys))
    monitor = 1

    Threads.@threads for i = 1:length(xys)
        zs[i] = loss_bomb(xys[i], bop)

        val = round(100*(monitor/length(xys)), sigdigits = 3)
        print(WHITE_BG("Surrogates: $(val)%   \r"))
        flush(stdout)
        monitor = monitor + 1
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

    optimizable = [bop.params[param].name for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]

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
        resolution = (2420, 2320),
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

    fig[-2,:] = Label(fig, L"[%$(bop.loss_func.name)] Loss(default) = %$(dl_ltx), Loss(opt) =  %$(ol_ltx)")

    fig[-1,:] = Label(fig, L"Defaults: $\delta =$ %$(δ_def), $a =$ %$(a_def), $\beta =$ %$(β_def), $A_\text{spring} =$ %$(A_spring_def), $\mu_\text{max} =$ %$(μ_max_def), $m =$ %$(m_def), $k_N =$ %$(k_N_def)")
    
    fig[0,:] = Label(fig, L"Optimals: $\delta =$ %$(δ_opt), $a =$ %$(a_opt), $\beta =$ %$(β_opt), $A_\text{spring} =$ %$(A_spring_opt), $\mu_\text{max} =$ %$(μ_max_opt), $m =$ %$(m_opt), $k_N =$ %$(k_N_opt)")

    outfile = joinpath(@__DIR__, "..", "figures", "opt_test.png")
    rm(outfile, force = true)
    save(outfile, fig)

    return fig
end