"""
    const OPTIMIZATION_PARAMETER_NAMES

A `Vector` of `String`s giving the names of all the parameters it is possible to optimize by default.

Equal to `["δ", "τ", "σ", "A_spring", "λ", "μ_max", "m", "k_N"]`.

### Definitions

- `δ`: The buoyancy of a clump. 
- `τ`: Measures the inertial response time of the medium to the particle.
- `σ`: A prefactor multiplying the Stokes drift term. 
- `A_spring`: The maximum value of the spring stiffness function.
- `λ`: A prefactor multiplying the spring natural length.
- `μ_max`:  Sargassum maximum growth rate.
- `m`: Sargassum mortality rate.
- `k_N`: Sargassum nutrient (N) uptake half saturation.
"""
const OPTIMIZATION_PARAMETER_NAMES = ["δ", "τ", "σ", "A_spring", "λ", "μ_max", "m", "k_N"]


"""
    struct LossFunction

A container for the function used for measuring the loss of a simulation.

### Fields 

- `f`: A `Function`. This function must be callable as `f(rtr::RaftTrajectory)` and return a `Real`. This measures the 
    "score" of a `RaftTrajectory` (lower is better).
- `metric`: A `Function`. This function must be callable as `f(a::Matrix, b::Matrix)` and return a `Real`. This measures
    the "distance" between two matrices, e.g. `(a, b) -> sum(abs.(a - b))` is the `L1` norm.
- `name`: A `String` giving the name of the loss function, e.g. `"L1"`.

### Constructor

`LossFunction(ymw1, ymw2, dists; corners, metric, cumulative, name)`

### Arguments 

- `ymw1`: A `NTuple{3, Integer}` giving the initial `(year, month, week)` of the integration.
- `ymw2`: A `NTuple{3, Integer}` giving the final `(year, month, week)` of the integration.
- `dists`: A `Dict` mapping `(year, month)`s to `SargassumDistribution`s spanning the integration time period.

### Optional Arguments

- `corners`: The `(lon_min, lon_max, lat_min, lat_max)` coordinates to restrict the loss calculation to. Default `(-180, 180, -90, 90)`.
- `metric`: As defined above. Default `(a, b) -> sum(abs.(a - b))`.
- `cumulative`: If `true`, the loss is added at each intermediate week of the integration rather than just at the end. Default `false`.
- `name`: A `String` giving the name of the loss function. Default `"L1"`.
"""
struct LossFunction
    f::Function
    metric::Function
    name::String

    function LossFunction(
        ymw1::NTuple{3, Integer}, 
        ymw2::NTuple{3, Integer},
        dists::Dict{Tuple{Int64, Int64}, SargassumFromAFAI.SargassumDistribution};
        corners::Tuple{Real, Real, Real, Real} = (-180, 180, -90, 90),
        metric::Function = (a, b) -> sum(abs.(a - b)), 
        cumulative::Bool = false,
        name::String = "L1")

        ymws = ymwspan2weekspan(ymw1, ymw2)
        tspans = [(ymw2time(ymws[i]...), ymw2time(ymws[i + 1]...)) for i = 1:length(ymws) - 1]
        # integrate for a week, then compare with the distribution (i.e. the next week's distribution)

        function weekly_loss(rtr::RaftTrajectory)
            loss_total = 0.0

            i1 = cumulative ? 1 : length(ymws)-1

            for i = i1:length(ymws)-1
                year, month, week = ymws[i+1]
                target = dists[(year, month)]

                lons = findall(lon -> corners[1] <= lon <= corners[2], target.lon)
                lats = findall(lat -> corners[3] <= lat <= corners[4], target.lat)

                data = bins(time_slice(rtr, tspans[i]), target) |> x -> x/sum(x)
                target = target.sargassum[:,:,week] |> x -> x/sum(x)
                loss_total = loss_total + metric(target[lons, lats], data[lons, lats])
            end

            return loss_total
        end


        return new(weekly_loss, metric, name)
    end
end


"""
    mutable struct OptimizationParameter{T}

A container for the data, values and bounds for a parameter to be optimized.

### Fields 

- `name`: A `String` with the name of the parameter, must be one of [`OPTIMIZATION_PARAMETER_NAMES`](@ref).
- `default`: The default (starting) value of the parameter.
- `bounds`: A `Tuple` giving the upper and lower bounds of the parameter.
- `opt`: A optimal value of the parameter.
- `optimizable`: A `Bool` such that, if `false`, the parameter will not be optimized away from its default value.

### Constructor

Use 

`OptimizationParameter(name, default, bounds, optimizable; opt = nothing)`
"""
mutable struct OptimizationParameter{T<:Real} 
    name::String
    default::T
    bounds::Tuple{T, T}
    opt::Union{Nothing, T}
    optimizable::Bool

    function OptimizationParameter(
        name::String, 
        default::Real, 
        bounds::Tuple{Real, Real},
        optimizable::Bool;
        opt::Union{Nothing, Real} = nothing)

        @assert name in OPTIMIZATION_PARAMETER_NAMES "Got $(name) ∉ $(OPTIMIZATION_PARAMETER_NAMES)"
        @assert first(bounds) < last(bounds)

        def, lb, ub = promote(default, first(bounds), last(bounds))

        op_opt = opt === nothing ? nothing : promote(opt, def)[1]

        return new{typeof(def)}(name, def, (lb, ub), op_opt, optimizable)
    end
end

"""
    mutable struct BOMBOptimizationProblem{T, U, C}

A container for all the data defining an optimization problem.

### Fields 

- `params`: A `Dict` mapping each element of `[OPTIMIZATION_PARAMETER_NAMES](@ref)` to an [`OptimizationParameter`](@ref) \
that contains it.
- `rhs`: The `Function` to integrate, generally should be [`Raft!`](@ref), but [`Leeway!`] can be used \
for testing purposes.
- `immortal`: A `Bool` such that if `true`, the [`ImmortalModel`](@ref) will be used, resulting in no clump \
growths or deaths.
- `ics`: The [`InitialConditions`](@ref) for the integration.
- `springs`: An [`BOMBSpring`](@ref) for the integration.
- `connections`: A subtype of [`AbstractConnections`](@ref) for the integration.
- `loss_func`: The [`LossFunction`](@ref) used during the optimization.
- `opt`: The minimal [`LossFunction`](@ref) obtained. If `nothing`, the problem is considered unoptimized.
- `opt_rtr`: When optimized, the [`RaftTrajectory`](@ref) that attained the optimal loss.
- `seed`: A `Random.seed!` used in the integration.

### Constructor

`BOMBOptimizationProblem(; params, rhs, immortal, ics, springs, connections, loss_func, seed = 1234)`

The fields `opt` and `opt_rtr` are initialized to `nothing`.
"""
mutable struct BOMBOptimizationProblem{T<:Real, U<:Integer, C<:AbstractConnections}
    params::Dict{String, OptimizationParameter{T}}
    rhs::Function
    immortal::Bool
    ics::InitialConditions{T}
    springs::BOMBSpring{T}
    connections::C
    loss_func::LossFunction
    opt::Union{Nothing, T}
    opt_rtr::Union{Nothing, RaftTrajectory{U, T}}
    seed::U

    function BOMBOptimizationProblem(;
        params::Dict{String, OptimizationParameter{T}},
        rhs::Function,
        immortal::Bool,
        ics::InitialConditions{T},
        springs::BOMBSpring{T},
        connections::C,
        loss_func::LossFunction,
        seed::U = 1234) where {T<:Real, U<:Integer, C<:AbstractConnections}

        if length(params) > 0 
            @warn "No parameters are set to be optimized."
        end
        @assert rhs in [Raft!, Leeway!]

        return new{T, U, C}(params, rhs, immortal, ics, springs, connections, loss_func, nothing, nothing, seed)
    end
end

"""
    optimizable(bop)

Return a vector of names of optimizable parameters of `bop`. These are ordered as they appear in `OPTIMIZATION_PARAMETER_NAMES`.
"""
function optimizable(bop::BOMBOptimizationProblem)
    return [bop.params[param].name for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]
end

"""
    RaftParameters(bop, type)

Construct [`RaftParameters`](@ref) for `bop::BOMBOptimizationProblem` according to `type`.

### Type

Type can be a `String` or a `Vector`. If it is a string, it can take one of two values:

- `"default"`: Each parameter is set equal to its `default`.
- `"opt"`: Each parameter is set equal to its optimal value, if it is both optimizable and optimized. Otherwise, `default` is used.

If it is a vector, it must have a length equal to `length(optimizable(bop))`. The `RaftParameters` is constructed 
by setting the optimizable parameters equal to the values in the vector (ordered as in `OPTIMIZATION_PARAMETER_NAMES`)
and leaving all other parameters equal to the defaults.
"""
function RaftParameters(bop::BOMBOptimizationProblem, type::Union{String, Vector{<:Real}})
    seed!(bop.seed)

    ps = OPTIMIZATION_PARAMETER_NAMES
    if type isa String
        @assert type in ["default", "opt"]   

        if type == "default"
            δ, τ, σ, A_spring, λ, μ_max, m, k_N = [bop.params[param].default for param in ps]
        elseif type == "opt"
            if bop.opt === nothing
                @warn "The `opt` parameter value has been selected, but `bop` has not been optimized. The default parameters will be used."
                δ, τ, σ, A_spring, λ, μ_max, m, k_N = [bop.params[param].default for param in ps]
            else
                δ, τ, σ, A_spring, λ, μ_max, m, k_N = [bop.params[param].optimizable ? bop.params[param].opt : bop.params[param].default for param in ps]
            end
        end
    elseif type isa Vector
        @assert length(type) == length(optimizable(bop)) "The vector `type` must have the same number of entries as the number of optimizable parameters."

        ps = OPTIMIZATION_PARAMETER_NAMES
        vs = []
        j = 1
        for i = 1:length(ps)
            if bop.params[ps[i]].optimizable
                push!(vs, type[j])
                j = j + 1
            else
                push!(vs, bop.params[ps[i]].default)
            end
        end

        δ, τ, σ, A_spring, λ, μ_max, m, k_N = vs
    end
    
    # ics
    ics = bop.ics |> deepcopy

    # clumps
    clumps = ClumpParameters(δ = δ, σ = σ)
    clumps = ClumpParameters(clumps.α, τ, clumps.R, clumps.f, clumps.σ)
    
    # springs
    springs = BOMBSpring(A_spring, λ*bop.springs.L)
    
    # connections
    connections = bop.connections |> deepcopy
    
    # growth-death
    if bop.immortal
        gd_model = ImmortalModel()
    else
        bmp = BrooksModelParameters(TEMPERATURE_ITP.x, NUTRIENTS_ITP.x, 
            clumps_limits = (0, Integer(2*ics.ics[1])), # the number of clumps can at most double
            μ_max = μ_max,
            m = m,
            k_N = k_N)
        gd_model = BrooksModel(params = bmp)
    end

    # land
    land = Land()
    
    return RaftParameters(
        ics = ics,
        clumps = clumps,
        springs = springs,
        connections = connections,
        gd_model = gd_model,
        land = land)
end

"""
    simulate(bop::BOMBOptimizationProblem, type; high_accuracy, showprogress)

Integrate `bop` by constructing the [`RaftParameters`](@ref) implied by its fields and calling [`simulate(::RaftParameters)`](@ref), 
returning a [`RaftTrajectory`](@ref).

### Type

Type can be a `String` or a `Vector`. If it is a string, it can take one of two values:

- `"default"`: Each parameter is set equal to its `default`.
- `"opt"`: Each parameter is set equal to its optimal value, if it is both optimizable and optimized. Otherwise, `default` is used.

If it is a vector, it must have a length equal to `length(optimizable(bop))`. The `RaftParameters` is constructed 
by setting the optimizable parameters equal to the values in the vector (ordered as in `OPTIMIZATION_PARAMETER_NAMES`)
and leaving all other parameters equal to the defaults.

### Optional Arguments

- `ics`: If provided, `bop` will be integrated with these `InitialConditions` rather than `bop.ics`. This has 
the effect of create a new `BOMBOptimizationProblem` with different `ics` and `springs`. The springs 
are updated with `L = ΔL(ics)`.
- `high_accuracy`: A `Bool` which, if `true`, uses higher tolerances in the integration. Default `false`.
- `showprogress`: A `Bool` which outputs the integration progress when `true`. Default `false`.
"""
function simulate(
    bop::BOMBOptimizationProblem,
    type::Union{String, Vector{<:Real}};
    ics::Union{InitialConditions, Nothing} = nothing,
    high_accuracy::Bool = false, 
    showprogress::Bool = false)

    seed!(bop.seed)

    if ics === nothing
        rp = RaftParameters(bop, type)
    else
        if "A_spring" in optimizable(bop) && type == "opt"
            A_spring = bop.params["A_spring"].opt
        else
            A_spring = bop.params["A_spring"].default
        end

        springs = BOMBSpring(A_spring, ΔL(ics)) # don't apply λ, might be nonsensical in general
        bop_new = deepcopy(bop)
        bop_new.ics = ics
        bop_new.springs = springs

        rp = RaftParameters(bop_new, type)
    end

    abstol = high_accuracy ? 1.0e-6 : nothing
    reltol = high_accuracy ? 1.0e-6 : nothing

    return simulate(rp, rhs = bop.rhs, abstol = abstol, reltol = reltol, showprogress = showprogress)
end


"""
    optimize!(bop; time_limit, high_accuracy, verbose)

Optimize the [`BOMBOptimizationProblem`](@ref) in `bop` using the `Metaheuristics` optimization package and update it
to include the optimal results. The Evolutionary Centers Algorithm is used via `Metaheuristics.ECA`.

### Arguments 

- `bop`:: An unoptimized [`BOMBOptimizationProblem`](@ref).

### Optional Arguments 

- `time_limit`: A `Float64` giving the upper time limit in seconds on the length of the optimization. Default `300.0`.
- `high_accuracy`: A `Bool` which, if `true`, uses higher tolerances in the integration. Default `false`.
- `verbose`: Show simplified results each iteration of the optimization. Default `true`.
"""
function optimize!(
    bop::BOMBOptimizationProblem; 
    time_limit::Float64 = 300.0,
    high_accuracy::Bool = false,
    verbose = true)

    seed!(bop.seed)
    @assert bop.opt === nothing "`bop` already has an optimal value"

    lb = [bop.params[param].bounds[1] for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]
    ub = [bop.params[param].bounds[2] for param in OPTIMIZATION_PARAMETER_NAMES if bop.params[param].optimizable]
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)

    function f_parallel(X)
        fitness = zeros(size(X,1))
        Threads.@threads for i in 1:size(X,1)
            fitness[i] = bop.loss_func.f(simulate(bop, X[i,:], high_accuracy = high_accuracy))
        end
        return fitness
    end

    options = Metaheuristics.Options(time_limit = time_limit, parallel_evaluation = true, verbose = verbose)
    algorithm = Metaheuristics.ECA(options = options)
    result = Metaheuristics.optimize(f_parallel, bounds, algorithm)

    bop.opt = result.best_sol.f
    i = 1
    for param in optimizable(bop)
        bop.params[param].opt = result.best_sol.x[i]
        i = i + 1
    end

    bop.opt_rtr = simulate(bop, "opt", high_accuracy = high_accuracy, showprogress = false)

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
- `high_accuracy`: A `Bool` which, if `true`, uses lower tolerances in the integration. Default `false`.
- `verbose`: Show the progress of the sampling. Default `true`.
"""
function sample!(
    bop::BOMBOptimizationProblem,
    n_samples::Integer;
    sampling_algorithm::QuasiMonteCarlo.SamplingAlgorithm = SobolSample(),
    time_limit::Float64 = 300.0,
    high_accuracy::Bool = false,
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

        cur_loss = bop.loss_func.f(simulate(bop, samps[:,i], high_accuracy = high_accuracy))
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

    bop.opt_rtr = simulate(bop, "opt", high_accuracy = high_accuracy, showprogress = false)

    return nothing 
end