using DifferentialEquations

function f(du, u, p, t)
    for i in 1:length(u)
        du[i] = 0.3 * u[i]
    end
end

function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    1 - maximum(u)
end

function affect!(integrator)
    # u = integrator.u
    # resize!(integrator, length(u) + 1)
    # maxidx = findmax(u)[2]
    # Θ = rand()
    # u[maxidx] = Θ
    # u[end] = 1 - Θ

    println(length(integrator.u))
    deleteat!(integrator, 2)
    println(length(integrator.u))
    nothing
end


callback = ContinuousCallback(condition, affect!)
u0 = [0.2, 0.05, 0.05, 0.05]
tspan = (0.0, 10.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), callback = callback)

sol.u .|> length