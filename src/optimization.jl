"""
    struct LossFunction{F1, F2}

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
struct LossFunction{F1<:Function, F2<:Function}
    f::F1
    metric::F2
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

            return isnan(loss_total) ? Inf : loss_total
        end


        return new{typeof(weekly_loss), typeof(metric)}(weekly_loss, metric, name)
    end
end


"""
    optimize!(f, param_bounds, loss_func; time_limit, verbose, seed)

Find the minimum of `f(X)` where `X ∈ param_bounds`. The Evolutionary Centers Algorithm is used via `Metaheuristics.ECA`.

### Arguments 

- `f`: The function to optimize. Should be evaluable at a vector whose length is equal to `length(param_bounds)` and return a `Real`.
- `param_bounds`: A vector of tuples giving the box constraints of the optimization.

### Optional Arguments 

- `time_limit`: A `Float64` giving the upper time limit in seconds on the length of the optimization. Default `300.0`.
- `verbose`: Show simplified results each iteration of the optimization. Default `true`.
- `seed`: An integer seed for randomness. Default `1234`.
"""
function optimize!(
    f::Function,
    param_bounds::Vector{<:Tuple{Real, Real}}; 
    time_limit::Float64 = 300.0,
    verbose = true,
    seed::Integer = 1234)

    seed!(seed)

    lb = [p[1] for p in param_bounds]
    ub = [p[2] for p in param_bounds]
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)

    function f_parallel(X)
        fitness = zeros(size(X,1))
        Threads.@threads for i in 1:size(X,1)
            fitness[i] = f(X[i,:])
        end
        return fitness
    end

    options = Metaheuristics.Options(time_limit = time_limit, parallel_evaluation = true, verbose = verbose)
    algorithm = Metaheuristics.ECA(options = options)
    result = Metaheuristics.optimize(f_parallel, bounds, algorithm)

    return result
end

# """
#     sample!(bop, n_samples; sampling_algorithm, time_limit, high_accuracy, verbose)

# Optimize the [`BOMBOptimizationProblem`](@ref) in `bop` using the `QuasiMonteCarlo` package. Specifically, 
# take `n_samples` using `QuasiMonteCarlo.SobolSample()` over all optimizable parameters so that the optimal 
# value is just the sample with the lowest loss function.

# Return (samples, losses).

# ### Arguments 

# - `bop`: An unoptimized [`BOMBOptimizationProblem`](@ref).
# - `n_samples`: The number of samples to take.

# ### Optional Arguments 

# - `sampling_algorithm`: A `QuasiMonteCarlo.SamplingAlgorithm`. Default `SobolSample()`.
# - `time_limit`: If provided, a `Float64` giving the upper time limit in seconds on the length of the optimization. Default `nothing`.
# - `high_accuracy`: A `Bool` which, if `true`, uses lower tolerances in the integration. Default `false`.
# - `verbose`: Show the progress of the sampling. Default `true`.
# """
# function sample!(
#     bop::BOMBOptimizationProblem,
#     n_samples::Integer;
#     sampling_algorithm::QuasiMonteCarlo.SamplingAlgorithm = SobolSample(),
#     time_limit::Union{Float64, Nothing} = nothing,
#     high_accuracy::Bool = false,
#     verbose = true)

#     lb = [bop.params[p].bounds[1] for p in OPTIMIZATION_PARAMETER_NAMES if bop.params[p].optimizable]
#     ub = [bop.params[p].bounds[2] for p in OPTIMIZATION_PARAMETER_NAMES if bop.params[p].optimizable]
#     samps = QuasiMonteCarlo.sample(n_samples, lb, ub, sampling_algorithm)

#     fitness = zeros(n_samples)
#     cur_best = Inf
#     timedout = false
#     tstart = time()

#     j = 0
#     Threads.@threads for i in 1:n_samples
#         j = j + 1
#         if (time_limit !== nothing) && (time() - tstart > time_limit)
#             timedout = true
#             break
#         end

#         cur_loss = bop.loss_func.f(simulate(bop, samps[:,i], high_accuracy = high_accuracy))
#         fitness[i] = cur_loss

#         if cur_loss < cur_best
#             cur_best = cur_loss
#         end

#         if verbose
#             val = 100*j/n_samples
#             print(WHITE_BG("Sampling: $(val)%, Best: $(cur_best) \r"))
#             flush(stdout)
#         end
#     end

#     if timedout
#         @info "Exceeded time limit of $(round(time_limit, sigdigits = 3)) seconds."
#     end

#     idx = findall(x -> x != 0.0, fitness)
#     fitness = fitness[idx]
#     samps = samps[:,idx]

#     best_loss, idx = findmin(fitness)
#     best_params = samps[:,idx]

#     bop.opt = best_loss
#     i = 1
#     for param in [name for name in OPTIMIZATION_PARAMETER_NAMES if bop.params[name].optimizable]
#         bop.params[param].opt = best_params[i]
#         i = i + 1
#     end

#     bop.opt_rtr = simulate(bop, "opt", high_accuracy = high_accuracy, showprogress = false)

#     return (samps, fitness)
# end