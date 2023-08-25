using DifferentialEquations

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 10.0)
prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(prob, Tsit5())

function f_extra(u)
    return sum(u)
end

function lorenz_plus!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    du[4] = f_extra(u[1:3])
end

u0 = [1.0; 0.0; 0.0]
u0_plus = [u0 ; f_extra(u0)]
tspan_plus = (0.0, 10.0)
prob_plus = ODEProblem(lorenz_plus!, u0_plus, tspan_plus)
sol_plus = solve(prob_plus, Tsit5())

################################################################

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 10.0)
prob = ODEProblem(lorenz!, u0, tspan)

integ = init(prob, Tsit5())


mutable struct ODEHolder{T<:Real}
    u::Vector{T}
end

function f_de!(du, u, p, t)
    du[1] = sum(p.u)
end

u_extra = [f_extra(u0)]
holder = ODEHolder(Float64[])
prob_extra = ODEProblem(f_de!, u_extra, tspan, holder)

integ_plus = init(prob_extra, Tsit5())



