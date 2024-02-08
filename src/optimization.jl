"""
    const OPTIMIZATION_PARAMETER_NAMES

A `Vector` of `String`s giving the names of all the parameters it is possible to optimize by default.

Equal to `["δ", "a", "σ", "A_spring", "λ", "μ_max", "m", "k_N"]`.
"""
const OPTIMIZATION_PARAMETER_NAMES = ["δ", "a", "σ", "A_spring", "λ", "μ_max", "m", "k_N"]

"""
    struct LossFunction

A container for a function used for measuring the distance between two binned histograms for \
optimization. 

### Fields 

- `f`: A `Function`. This function must be callable as `f(a, b)`, where `a` and `b` are matrices and it \
must return a `Real` such that `f(a, a) = 0.0`. Example: `(a, b) -> sum((a - b) .^ 2)`. If the loss function 
needs to depend on the exact form of `a` and `b`, it is assumed that `a` is the target distribution and `b` is the
distribution being tested. 
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
const LOSS_L1 = LossFunction(f = (a::Matrix, b::Matrix) -> sum(abs.(a - b)) , name = "L1")

"""
    const LOSS_L2

A `LossFunction` for the L2 norm.
"""
const LOSS_L2 = LossFunction(f = (a::Matrix, b::Matrix) -> sqrt(sum((a - b) .^ 2)) , name = "L2")

"""
    const LOSS_COR

A `LossFunction` for the negative correlation.
"""
const LOSS_COR = LossFunction(f = (a::Matrix, b::Matrix) -> -cor(vec(a), vec(b)), name = "-COR")

function l2_shape(a::Matrix, b::Matrix)
    # a is target, b is data
    idx = findall(iszero, a) 
    va, vb = a[idx], b[idx]
    return sqrt(sum((va - vb) .^ 2))
end

"""
    const LOSS_L2_SHAPE

A `LossFunction` for the L2 norm, calculated only at those locations where the target is zero. In other words, 
the loss is minimized when the data does not show up where it shouldn't (but doesn't necessarily show up where it should).
"""
const LOSS_L2_SHAPE = LossFunction(f = l2_shape, name = "L2_SHAPE")

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
- `opt_rtr`: When optimized, the [`RaftTrajectory`](@ref) that attained the optimal loss.
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
    opt_rtr::Union{Nothing, RaftTrajectory{U, T}}
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
        opt_rtr::Union{Nothing, RaftTrajectory{U, T}} = nothing,
        seed::U = 1234) where {T<:Real, U<:Integer}

        @assert length(params) > 0 "Must optimize at least one parameter."
        @assert rhs in [Raft!, WaterWind!]
        @assert DateTime(first(tspan)...) < DateTime(last(tspan)...)
        @assert n_levels > 0 "Need at least one level"

        return new{T, U}(params, rhs, immortal, tspan, n_levels, t_extra, loss_func, opt, opt_rtr, seed)
    end
end

function Base.show(io::IO, bop::BOMBOptimizationProblem)
    optimizable = [bop.params[param].name for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]
    space = [(name, bop.params[name].default, (bop.params[name].bounds)) for name in optimizable]
    optimized_q = bop.opt === nothing ? false : true
    
    println(io, "BOMBOptimizationProblem")

    println(io, " Optimized?: $(optimized_q)")
    println(io, " Search space: $(space)")

    if optimized_q
        opts = Tuple{String, Float64}[]
        for name in OPTIMIZATION_PARAMETER_NAMES
            opt = bop.params[name].opt
            if opt === nothing 
                opt = bop.params[name].default
            end
            push!(opts, (name, round(Float64(opt), sigdigits = 4)))
        end
        println(io, " Optimals: $(opts)")
        println(io, " Loss: $(bop.opt)")
    end
end

"""
    simulate(bop::BOMBOptimizationProblem; high_accuracy, type, showprogress)

Integrate `bop` by constructing the [`RaftParameters`](@ref) implied by its fields using [`simulate(::RaftParameters)`](@ref).

### Optional Arguments

- `high_accuracy`: A `Bool` which, if `true`, uses higher tolerances in the integration. Default `true`.
- `type`: A `String` identifying the [`OptimizationParameter`](@ref) values to be used during the simulation.
    - `"val"`: The default case, each parameter is equal to its `val`; used during optimization.
    - `"default"`: Each parameter is set equal to its `default`.
    - `"opt`: Each parameter is set equal to its optimal value, if it is both optimizable and optimized. Otherwise, `default` is used.
- `showprogress`: A `Bool` which outputs the integration progress when `true`. Default `false`.
"""
function simulate(
    bop::BOMBOptimizationProblem;
    high_accuracy::Bool = true, 
    type::String = "val",
    showprogress::Bool = false)

    seed!(bop.seed)

    @assert type in ["val", "default", "opt"]

    if type == "val"
        δ, a, σ, A_spring, λ, μ_max, m, k_N = [bop.params[param].val for param in OPTIMIZATION_PARAMETER_NAMES]
    elseif type == "default"
        δ, a, σ, A_spring, λ, μ_max, m, k_N = [bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES]
    elseif type == "opt"
        if bop.opt === nothing
            @warn "The `opt` parameter value has been selected, but `bop` has not been optimized. The default parameters will be used."
            δ, a, σ, A_spring, λ, μ_max, m, k_N = [bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES]
        else
            δ, a, σ, A_spring, λ, μ_max, m, k_N = [bop.params[param].optimizable ? bop.params[param].opt : bop.params[param].default for param in OPTIMIZATION_PARAMETER_NAMES]
        end
    end

    # time
    start_date, end_date = bop.tspan
    dist = SargassumFromAFAI.DIST_2018[start_date]
    tspan = yearmonth2tspan(start_date, end_date, t_extra = (0, bop.t_extra))
    
    # clumps
    clumps = ClumpParameters(δ = δ, a = a, σ = σ)
    
    # springs
    L_spring = λ*ΔL(dist)
    function spring_k(x::Real; A::Real = A_spring, L::Real = L_spring)
        return A * (exp((x - 2*L)/0.2) + 1)^(-1)
    end

    springs = SpringParameters(spring_k, L_spring)
    
    # initial conditions
    ics = InitialConditions(dist, [1], bop.n_levels, "levels", EQR_DEFAULT)
    
    # connections
    connections = ConnectionsNearest(10)
    
    # growth-death
    if bop.immortal
        gd_model = ImmortalModel()
    else
        bmp = BrooksModelParameters(TEMP_ITP.x, NO3_ITP.x, 
            clumps_limits = (0, Integer(2*ics.ics[1])), # the number of clumps can at most double
            μ_max = μ_max,
            m = m,
            k_N = k_N)
        gd_model = BrooksModel(params = bmp)
    end

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
    
    abstol = high_accuracy ? 1.0e-6 : nothing
    reltol = high_accuracy ? 1.0e-6 : nothing

    return simulate(rp, rhs = bop.rhs, abstol = abstol, reltol = reltol, showprogress = showprogress)
end


"""
    loss(u::Vector{<:Real}, bop::BOMBOptimizationProblem; high_accuracy, type, showprogress)

Compute the loss associated with parameter values `u` in `bop`.

The vector `u` should be the same length as the number of `bop.params` which are optimizable. Then, 
the `bop` is integrated with the values of each parameter set equal to the entries of `u` in the 
order defined by [`OPTIMIZATION_PARAMETER_NAMES`](@ref).

### Optional Arguments

The following arguments are passed directly to `simulate`.

- `high_accuracy`: A `Bool` which, if `true`, uses higher tolerances in the integration. Default `true`.
- `type`: A `String` identifying the [`OptimizationParameter`](@ref) values to be used during the simulation.
    - `"val"`: The default case, each parameter is equal to its `val`; used during optimization.
    - `"default"`: Each parameter is set equal to its `default`.
    - `"opt`: Each parameter is set equal to its optimal value, if it is both optimizable and optimized. Otherwise, `default` is used.
- `showprogress`: A `Bool` which outputs the integration progress when `true`. Default `false`.
"""
function loss(
    u::Vector{<:Real}, 
    bop::BOMBOptimizationProblem;
    high_accuracy::Bool = true,
    type::String = "val",
    showprogress::Bool = false)
    optimizable = [bop.params[param].name for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]

    @assert length(u) == length(optimizable) "The vector u must have the same number of entries as the number of optimizable parameters."

    for i = 1:length(optimizable)
        bop.params[optimizable[i]].val = u[i]
    end

    start_date, end_date = bop.tspan
    tstart, tend = yearmonth2tspan(start_date, end_date, t_extra = (0, bop.t_extra))

    target_dist = SargassumFromAFAI.DIST_2018[end_date]
    target = target_dist.sargassum[:,:,1]
    target = target/sum(target)

    rtr = simulate(bop, high_accuracy = high_accuracy, type = type, showprogress = showprogress)
    rtr_final = time_slice(rtr, (tend - bop.t_extra, tend))
    data = bins(rtr_final, target_dist)
    data = data/sum(data)
    loss_val = bop.loss_func.f(target, data)

    if (type == "val") 
        if (bop.opt === nothing) || (loss_val < bop.opt)
            bop.opt = loss_val
            bop.opt_rtr = rtr
        end
    end
        
    return loss_val
end

"""
    loss(rtr::RaftTrajectory, bop::BOMBOptimizationProblem)

Compute the loss associated with the [`RaftTrajectory`](@ref) `rtr`.

Use this instead of `loss(u::Vector, bop)` when the integration is already done.

### Optional Arguments

None.
"""
function loss(
    rtr::RaftTrajectory, 
    bop::BOMBOptimizationProblem)

    start_date, end_date = bop.tspan
    tstart, tend = yearmonth2tspan(start_date, end_date, t_extra = (0, bop.t_extra))

    target_dist = SargassumFromAFAI.DIST_2018[end_date]
    target = target_dist.sargassum[:,:,1]
    target = target/sum(target)

    data = bins(time_slice(rtr, (tend - bop.t_extra, tend)), target_dist)
    data = data/sum(data)

    return bop.loss_func.f(target, data)
end

"""
    optimize!(bop; time_limit, high_accuracy, verbose)

Optimize the [`BOMBOptimizationProblem`](@ref) in `bop` using the `Metaheuristics` optimization package and update it
to include the optimal results. The Evolutionary Centers Algorithm is used via `Metaheuristics.ECA`.

### Arguments 

- `bop`:: An unoptimized [`BOMBOptimizationProblem`](@ref).

### Optional Arguments 

- `time_limit`: A `Float64` giving the upper time limit in seconds on the length of the optimization. Default `300.0`.
- `high_accuracy`: A `Bool` which, if `true`, uses lower tolerances in the integration. Default `true`.
- `verbose`: Show simplified results each iteration of the optimization. Default `true`.
"""
function optimize!(
    bop::BOMBOptimizationProblem; 
    time_limit::Float64 = 300.0,
    high_accuracy::Bool = true,
    verbose = true)

    @assert bop.opt === nothing "`bop` already has an optimal value"

    lb = [bop.params[param].bounds[1] for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]
    ub = [bop.params[param].bounds[2] for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)

    function f_parallel(X)
        fitness = zeros(size(X,1))
        Threads.@threads for i in 1:size(X,1)
            fitness[i] = loss(X[i,:], bop, high_accuracy = high_accuracy)
        end
        return fitness
    end

    options = Metaheuristics.Options(time_limit = time_limit, parallel_evaluation = true, verbose = verbose)
    algorithm = Metaheuristics.ECA(options = options)
    result = Metaheuristics.optimize(f_parallel, bounds, algorithm)

    bop.opt = result.best_sol.f
    i = 1
    for param in [name for name in OPTIMIZATION_PARAMETER_NAMES if bop.params[name].optimizable]
        bop.params[param].opt = result.best_sol.x[i]
        i = i + 1
    end

    return nothing
end

"""
    sample!(bop, n_samples; sampling_algorithm, time_limit, high_accuracy, verbose)

Optimize the [`BOMBOptimizationProblem`](@ref) in `bop` using the `QuasiMonteCarlo` package. Specifically, 
take `n_samples` using `QuasiMonteCarlo.SobolSample()` over all optimizable parameters so that the optimal 
value is just the sample with the lowest loss function.

### Arguments 

- `bop`: An unoptimized [`BOMBOptimizationProblem`](@ref).
- `n_samples`: The number of samples to take.

### Optional Arguments 

- `sampling_algorithm`: A `QuasiMonteCarlo.SamplingAlgorithm`. Default `SobolSample()`.
- `time_limit`: A `Float64` giving the upper time limit in seconds on the length of the optimization. Default `300.0`.
- `high_accuracy`: A `Bool` which, if `true`, uses lower tolerances in the integration. Default `true`.
- `verbose`: Show the progress of the sampling. Default `true`.
"""
function sample!(
    bop::BOMBOptimizationProblem,
    n_samples::Integer;
    sampling_algorithm::QuasiMonteCarlo.SamplingAlgorithm = SobolSample(),
    time_limit::Float64 = 300.0,
    high_accuracy::Bool = true,
    verbose = true)

    lb = [bop.params[p].bounds[1] for p in OPTIMIZATION_PARAMETER_NAMES if bop.params[p].optimizable]
    ub = [bop.params[p].bounds[2] for p in OPTIMIZATION_PARAMETER_NAMES if bop.params[p].optimizable]
    samps = QuasiMonteCarlo.sample(n_samples, lb, ub, sampling_algorithm)

    fitness = zeros(n_samples)
    cur_best = Inf
    timedout = false
    tstart = time()

    j = 0
    Threads.@threads for i in 1:n_samples
        j = j + 1
        if time() - tstart > time_limit
            timedout = true
            break
        end

        cur_loss = loss(samps[:,i], bop, high_accuracy = high_accuracy)
        fitness[i] = cur_loss

        if cur_loss < cur_best
            cur_best = cur_loss
        end

        if verbose
            val = 100*j/n_samples
            print(WHITE_BG("Sampling: $(val)%, Best: $(cur_best) \r"))
            flush(stdout)
        end
    end

    if timedout
        @info "Exceeded time limit of $(round(time_limit, sigdigits = 3)) seconds."
    end

    idx = findall(x -> x != 0.0, fitness)
    fitness = fitness[idx]
    samps = samps[:,idx]

    best_loss, idx = findmin(fitness)
    best_params = samps[:,idx]

    bop.opt = best_loss
    i = 1
    for param in [name for name in OPTIMIZATION_PARAMETER_NAMES if bop.params[name].optimizable]
        bop.params[param].opt = best_params[i]
        i = i + 1
    end

    return nothing 
end