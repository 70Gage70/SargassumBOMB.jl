include(joinpath(@__DIR__, "optimization.jl"))

# ["δ", "a", "β", "A_spring", "μ_max", "m", "k_N"]

###############################################################################################
# ## PHYSICS OPT

# ###
# initial_time = (2018, 3)
# final_time = (2018, 4)
# t_extra = 7

# δ_param = OptimizationParameter("δ",                1.25,   (1.05, 1.5),        true)
# a_param = OptimizationParameter("a",                1.0e-4, (1.0e-5, 1.0e-3),   true)
# β_param = OptimizationParameter("β",                0.0,    (0.0, 0.05),        true)
# A_spring_param = OptimizationParameter("A_spring",  1.0,    (0.1, 3.0),         true)
# μ_max_param = OptimizationParameter("μ_max",        0.1,    (0.05, 0.5),        false)
# m_param = OptimizationParameter("m",                0.05,   (0.0, 0.1),         false)
# k_N_param = OptimizationParameter("k_N",            0.012,  (0.005, 0.05),      false)

# ###
# bop = BOMBOptimizationProblem(
#     params = [δ_param, a_param, β_param, A_spring_param, μ_max_param, m_param, k_N_param],
#     rhs = Raft!,
#     immortal = true,
#     tspan = (initial_time, final_time),
#     n_levels = 10,
#     t_extra = 7,
# )

# ### 10/10: -0.5178
# n_sur_samples = 1000
# maxiters_opt = 10

# ###
# @info "Computing surrogate."
# @time rb = surrogate_bomb(n_sur_samples, bop)
# rm("data_physics_k_stokes_l1.jld2", force = true)
# jldsave("data_physics_k_stokes_l1.jld2"; bop = bop, rb = rb)

###############################################################################################
### BIOLOGY (with physics) OPT

# opt_data = jldopen(joinpath(@__DIR__, "data_physics_k.jld2"))
# rb_data = opt_data["rb"]
# sp = sortperm(rb_data.y)
# x_opt_phys = rb_data.x[sp][1] # δ, a, β, A_spring

# ###
# initial_time = (2018, 3)
# final_time = (2018, 4)
# t_extra = 7

# δ_param = OptimizationParameter("δ",                x_opt_phys[1],   (1.05, 1.5),        false)
# a_param = OptimizationParameter("a",                x_opt_phys[2],   (1.0e-5, 1.0e-3),   false)
# β_param = OptimizationParameter("β",                x_opt_phys[3],   (0.0, 0.01),        false)
# A_spring_param = OptimizationParameter("A_spring",  x_opt_phys[4],   (0.1, 3.0),         false)
# μ_max_param = OptimizationParameter("μ_max",        0.1,        (0.05, 0.5),        true)
# m_param = OptimizationParameter("m",                0.05,       (0.0, 0.1),         true)
# k_N_param = OptimizationParameter("k_N",            0.012,      (0.005, 0.05),      true)

# ###
# bop = BOMBOptimizationProblem(
#     params = [δ_param, a_param, β_param, A_spring_param, μ_max_param, m_param, k_N_param],
#     rhs = Raft!,
#     immortal = false,
#     tspan = (initial_time, final_time),
#     n_levels = 10,
#     t_extra = 7,
# )

# ###
# n_sur_samples = 200
# maxiters_opt = 10

# ###
# @info "Computing surrogate."
# @time rb = surrogate_bomb(n_sur_samples, bop)
# rm("data_bio_k.jld2", force = true)
# jldsave("data_bio_k.jld2"; bop = bop, rb = rb)
###
# @info "Optimizing surrogate."
# @time bop = optimize_bomb(rb, bop, maxiters_opt = maxiters_opt)
# @info "Plotting."
# @time plot_bop(bop)

###############################################################################################
## PHYSICS PLOT

opt_data = jldopen(joinpath(@__DIR__, "data_physics_k_stokes_l1.jld2"))
rb_data = opt_data["rb"]
sp = sortperm(rb_data.y)
x_opt_phys = rb_data.x[sp][1] # δ, a, β, A_spring

initial_time = (2018, 3)
final_time = (2018, 4)
t_extra = 7

δ_param = OptimizationParameter("δ",                1.25,   (1.05, 1.5),        true, opt = x_opt_phys[1])
a_param = OptimizationParameter("a",                1.0e-4, (1.0e-5, 1.0e-3),   true, opt = x_opt_phys[2])
β_param = OptimizationParameter("β",                0.0,    (0.0, 0.01),        true, opt = x_opt_phys[3])
A_spring_param = OptimizationParameter("A_spring",  1.0,    (0.1, 3.0),         true, opt = x_opt_phys[4])
μ_max_param = OptimizationParameter("μ_max",        0.1,    (0.05, 0.5),        false)
m_param = OptimizationParameter("m",                0.05,   (0.0, 0.1),         false)
k_N_param = OptimizationParameter("k_N",            0.012,  (0.005, 0.05),      false)

bop = BOMBOptimizationProblem(
    params = [δ_param, a_param, β_param, A_spring_param, μ_max_param, m_param, k_N_param],
    rhs = Raft!,
    immortal = true,
    tspan = (initial_time, final_time),
    n_levels = 10,
    t_extra = 7,
    opt = rb_data.y[sp][1]
)

@time plot_bop(bop)

##############################################################################################
### PHYSICS + BIO PLOT

# opt_data = jldopen(joinpath(@__DIR__, "data_physics_k.jld2"))
# rb_data = opt_data["rb"]
# sp = sortperm(rb_data.y)
# x_opt_phys = rb_data.x[sp][1] # δ, a, β, A_spring

# opt_data = jldopen(joinpath(@__DIR__, "data_bio_k.jld2"))
# rb_data = opt_data["rb"]
# sp = sortperm(rb_data.y)
# x_opt_bio = rb_data.x[sp][1] # μ_max, m, k_N

# initial_time = (2018, 3)
# final_time = (2018, 4)
# t_extra = 7

# δ_param = OptimizationParameter("δ",                1.25,   (1.05, 1.5),        true, opt = x_opt_phys[1])
# a_param = OptimizationParameter("a",                1.0e-4, (1.0e-5, 1.0e-3),   true, opt = x_opt_phys[2])
# β_param = OptimizationParameter("β",                0.0,    (0.0, 0.01),        true, opt = x_opt_phys[3])
# A_spring_param = OptimizationParameter("A_spring",  1.0,    (0.1, 3.0),         true, opt = x_opt_phys[4])
# μ_max_param = OptimizationParameter("μ_max",        0.1,    (0.05, 0.5),        true, opt = x_opt_bio[1])
# m_param = OptimizationParameter("m",                0.05,   (0.0, 0.1),         true, opt = x_opt_bio[2])
# k_N_param = OptimizationParameter("k_N",            0.012,  (0.005, 0.05),      true, opt = x_opt_bio[3])

# bop = BOMBOptimizationProblem(
#     params = [δ_param, a_param, β_param, A_spring_param, μ_max_param, m_param, k_N_param],
#     rhs = Raft!,
#     immortal = false,
#     tspan = (initial_time, final_time),
#     n_levels = 10,
#     t_extra = 7,
#     opt = rb_data.y[sp][1] # from data_bio.jld2
# )

# @time plot_bop(bop)