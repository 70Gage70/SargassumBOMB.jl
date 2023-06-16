using Interpolations, DifferentialEquations, ModelingToolkit
using LinearAlgebra: ⋅

include("parameters_mwe.jl")

lon = -180:1:180
lat = -90:0.5:90
time = 0.0:1:365
water = fill(1.0, length(lon), length(lat), length(time))
wind = fill(1.0, length(lon), length(lat), length(time))

interpolant_type = BSpline(Cubic(Line(OnGrid()))) 
const water_itp = scale(interpolate(water, interpolant_type), lon, lat, time)
const wind_itp = scale(interpolate(water, interpolant_type), lon, lat, time)

MaterialDerivativeX(VF, x, y, t) = gradient(VF, x, y, t) ⋅ [VF(x, y, t), VF(x, y, t), 1.0]
MaterialDerivativeY(VF, x, y, t) = gradient(VF, x, y, t) ⋅ [VF(x, y, t), VF(x, y, t), 1.0]
Vorticity(VF, x, y, t) = gradient(VF, x, y, t)[1] - gradient(VF, x, y, t)[2]

@variables t            
ddt = Differential(t)

v_x(x, y, t) = water_itp(x, y, t)
v_y(x, y, t) =  water_itp(x, y, t)
Dv_xDt(x, y, t) = MaterialDerivativeX(water_itp, x, y, t)
Dv_yDt(x, y, t) = MaterialDerivativeY(water_itp, x, y, t)
u_x(x, y, t, α) = (1 - α) * water_itp(x, y, t) + α * wind_itp(x, y, t)
u_y(x, y, t, α) = (1 - α) * water_itp(x, y, t) + α * wind_itp(x, y, t)
Du_xDt(x, y, t, α) = (1 - α) * MaterialDerivativeX(water_itp, x, y, t) + α * MaterialDerivativeX(wind_itp, x, y, t) 
Du_yDt(x, y, t, α) = (1 - α) * MaterialDerivativeY(water_itp, x, y, t) + α * MaterialDerivativeY(wind_itp, x, y, t) 
ω(x, y, t) = Vorticity(water_itp, x, y, t)

@register_symbolic v_x(x, y, t)
@register_symbolic v_y(x, y, t)
@register_symbolic Dv_xDt(x, y, t)
@register_symbolic Dv_yDt(x, y, t)
@register_symbolic u_x(x, y, t, α)
@register_symbolic u_y(x, y, t, α)
@register_symbolic Du_xDt(x, y, t, α)
@register_symbolic Du_yDt(x, y, t, α)
@register_symbolic ω(x, y, t)

@register_symbolic spring_force_x(x1, x2, y1, y2, parameters)
@register_symbolic spring_force_y(x1, x2, y1, y2, parameters)

function Clump(
    xy0::Vector{<:Real},
    clump_parameters::BOMParameters;
    name::Symbol)

    ps = @parameters α = clump_parameters.α τ = clump_parameters.τ R = clump_parameters.R f = clump_parameters.f
    @variables x(t) = xy0[1] y(t) = xy0[2]
    @variables Fx(t) Fy(t)

    eqs = [
    ddt(x) ~ u_x(x, y, t, α) + τ * (
        R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α) + Fx
    ),
    ddt(y) ~ u_y(x, y, t, α) + τ * (
        R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α) + Fy
    )
    ]

    return ODESystem(eqs, t, [x, y, Fx, Fy], ps; name)
end

function Raft(
    xy0::Matrix{<:Real},
    clump_parameters::Vector{<:BOMParameters},
    spring_parameters::Matrix{<:SpringParameters};
    name::Symbol)

    N_clumps = size(xy0, 1)
    @named clump 1:N_clumps i -> Clump(xy0[i,:], clump_parameters[i])

    forces_x = [clump[i].Fx ~ sum([spring_force_x(clump[i].x, clump[j].x, clump[i].y, clump[j].y, spring_parameters[i, j])] for j = 1:N_clumps if j != i)[1]
        for i = 1:N_clumps
    ]

    forces_y = [clump[i].Fy ~ sum([spring_force_y(clump[i].x, clump[j].x, clump[i].y, clump[j].y, spring_parameters[i, j])] for j = 1:N_clumps if j != i)[1]
        for i = 1:N_clumps
    ]

    eqs = [forces_x ; forces_y]

    return compose(ODESystem(eqs, t; name = name), clump[1:N_clumps]...)
end

function RectangularRaft(
    x_range::AbstractRange{<:Real}, 
    y_range::AbstractRange{<:Real},
    clump_parameters::BOMParameters, 
    spring_parameters::SpringParameters; 
    name::Symbol)

    # a rectangular mesh
    network = reverse.(collect(Iterators.product(reverse(y_range), x_range))) # reverse so that the first row has the largest y
    n_col = length(x_range)
    n_row = length(y_range)
    N_clumps = length(network)

    xy0 = Matrix{Float64}(undef, N_clumps, 2)
    clump_parameters_raft = Vector{BOMParameters}(undef, N_clumps)
    spring_parameters_raft = Matrix{SpringParameters}(undef, N_clumps, N_clumps)

    n(i, j) = (i - 1) * n_col + j

    for i = 1:n_row
        for j = 1:n_col
            xy0[n(i, j), :] .= network[i, j] 
            clump_parameters_raft[n(i, j)] = clump_parameters
            connections = filter(idx -> (1 <= idx[1] <= n_row) && (1 <= idx[2] <= n_col), [(i-1, j), (i+1, j), (i, j-1), (i, j+1)])
            connections = map(x->n(x...), connections)

            spring_parameters_raft[n(i, j), connections] .= spring_parameters
            spring_parameters_raft[n(i, j), setdiff(1:N_clumps, connections)] .= SpringParameters(k->0, 0)
        end
    end

    return Raft(xy0, clump_parameters_raft, spring_parameters_raft, name = name)
end

function sol(n)
    @info "Generating model."

    x0, y0 = (0.0, 0.0)
    x_range = range(start = x0 - 5, length = n, stop = x0 + 5)
    y_range = range(start = y0 - 5, length = n, stop = y0 + 5)
    clump_parameters = BOMParameters()
    spring_parameters = SpringParameters(k -> 20, 0.4)

    @named RRaft = RectangularRaft(x_range, y_range, clump_parameters, spring_parameters)
    @time RRaft = structural_simplify(RRaft)

    t_range = (0.0, 1.0)

    @info "Generating problem."

    @time prob = ODEProblem(RRaft, [], t_range, [], jac = true, sparse = true)

    @info "Solving problem."

    @time sol = solve(prob)

    return sol
end