using ProgressMeter
using MAT

include("models.jl")
include(joinpath(@__DIR__, "..", "plotting.jl"))

###########################################

function clump_t(n_traj_x, n_traj_y, t_short, t_long, t_span)
    xy0 = [0.0, 0.0] # arbitrary
    params = ClumpParameters(ref_itp)
    @named clump = Clump(xy0, params)
    clump = structural_simplify(clump)
    
    prob = ODEProblem(
        clump, 
        [],
        t_span, 
        []
    )
    
    # callbacks, so computations end when a particle exits the water
    condition(u, t, integrator) = (abs(water_itp.u(u..., t)) < 0.1) && (abs(water_itp.v(u..., t)) < 0.1)
    affect!(integrator) = terminate!(integrator)
    # save_positions ensures we don't save before or after the condition is satisfied, avoiding partial trajectories
    cb = DiscreteCallback(condition, affect!, save_positions = (false, false))  


    # generating initial conditions for each run
    # t_long = 100.0
    # t_short = 5.0
    # n_traj_x = 100
    # n_traj_y = 100

    n_steps = Integer(t_long/t_short) + 1
    x_range = range(start = -100.0 + 0.5, length = n_traj_x, stop = -50.0 - 0.5)
    y_range = range(start = 5.0 + 0.5, length = n_traj_y, stop = 35.0 - 0.5)
    x_range, y_range = sph2xy(x_range, y_range, ref_itp)
    xy0 = [
        [x, y] + [step(x_range)*(rand() - 1/2), step(y_range)*(rand() - 1/2)]  # random wiggle in x0, y0 to probe space more
        for x in x_range for y in y_range
    ] 

    n_traj = length(xy0)

    # prob_fun tells the ensemble how to "start over", in this case the i'th run is the problem remade with the i'th initial condition and time range.
    function prob_func(prob, i, repeat)
        remake(prob, u0 = xy0[i])
    end
    
    # output_func generally handles what is done with the solution at each iteration; here we only use it to keep track of progress
    progress = Progress(n_traj)
    
    function output_func(sol, i)
        next!(progress)
        sol, false
    end
    
    Eprob = EnsembleProblem(
        prob, 
        prob_func = prob_func,
        output_func = output_func
    )
    
    # @info "Solving model."
    
    sim = solve(Eprob, EnsembleThreads(), trajectories = n_traj, saveat = t_short, callback = cb);
    
    finish!(progress)    

    x0_long_t = Float64[]
    y0_long_t = Float64[]
    xT_long_t = Float64[]
    yT_long_t = Float64[]
    
    x0_short_t = Float64[]
    y0_short_t = Float64[]
    xT_short_t = Float64[]
    yT_short_t = Float64[]
    
    for sol in sim.u
        if length(sol.t) == n_steps
            xy_sol = xy2sph(sol.u, ref_itp)
    
            push!(x0_long_t, xy_sol[1, 1])
            push!(y0_long_t, xy_sol[1, 2])
            push!(xT_long_t, xy_sol[end, 1])
            push!(yT_long_t, xy_sol[end, 2])
            
            for j = 1:n_steps - 1
                push!(x0_short_t, xy_sol[j, 1])
                push!(y0_short_t, xy_sol[j, 2])
                push!(xT_short_t, xy_sol[j + 1, 1])
                push!(yT_short_t, xy_sol[j + 1, 2])
            end
        end
    end    

    return (x0_long_t, y0_long_t, xT_long_t, yT_long_t, x0_short_t, y0_short_t, xT_short_t, yT_short_t)
end

function all_clumps(n_traj_x, n_traj_y, t_short, t_long, t_stop)
    x0_long = Float64[]
    y0_long = Float64[]
    xT_long = Float64[]
    yT_long = Float64[]

    x0_short = Float64[]
    y0_short = Float64[]
    xT_short = Float64[]
    yT_short = Float64[]

    t_spans = [(i, i + 100.0) for i = 0.0:1.0:t_stop]

    for i = 1:length(t_spans)
        x0_long_t, y0_long_t, xT_long_t, yT_long_t, x0_short_t, y0_short_t, xT_short_t, yT_short_t = clump_t(n_traj_x, n_traj_y, t_short, t_long, t_spans[i])
        
        x0_long = [x0_long ; x0_long_t]
        y0_long = [y0_long ; y0_long_t]
        xT_long = [xT_long ; xT_long_t]
        yT_long = [yT_long ; yT_long_t]
    
        x0_short = [x0_short ; x0_short_t]
        y0_short = [y0_short ; y0_short_t]
        xT_short = [xT_short ; xT_short_t]
        yT_short = [yT_short ; yT_short_t]
        
        println(i/length(t_spans))
    end

    matwrite("x0xT_BOM_clump_5days.mat", Dict("x0" => x0_short, "y0" => y0_short, "xT" => xT_short, "yT" => yT_short))
    matwrite("x0xT_BOM_clump_100days.mat", Dict("x0" => x0_long, "y0" => y0_long, "xT" => xT_long, "yT" => yT_long))

    return nothing
end

        
