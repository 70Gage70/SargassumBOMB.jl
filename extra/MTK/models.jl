using DifferentialEquations
using ModelingToolkit
using JLD2

include(joinpath(@__DIR__, "..", "parameters.jl"))
include(joinpath(@__DIR__, "..", "..", "interpolants", "itp-derivatives.jl"))

# itp_path = joinpath(@__DIR__, "..", "..", "interpolants", "ias")
# isdefined(@__MODULE__, :water_itp) || (const water_itp = load(joinpath(itp_path, "water_itp.jld2"), "water_itp"))
# isdefined(@__MODULE__, :wind_itp) || (const wind_itp = load(joinpath(itp_path, "wind_itp.jld2"), "wind_itp"))
# isdefined(@__MODULE__, :ref_itp) || (const ref_itp = water_itp.ref)

itp_path = joinpath(@__DIR__, "..", "..", "interpolants", "glorys")
isdefined(@__MODULE__, :water_itp) || (const water_itp = load(joinpath(itp_path, "water_itp.jld2"), "water_itp"))
isdefined(@__MODULE__, :wind_itp) || (const wind_itp = load(joinpath(itp_path, "wind_itp.jld2"), "wind_itp"))
isdefined(@__MODULE__, :temp_itp) || (const temp_itp = load(joinpath(itp_path, "temp_itp.jld2"), "temp_itp"))
isdefined(@__MODULE__, :no3_itp) || (const no3_itp = load(joinpath(itp_path, "no3_itp.jld2"), "no3_itp"))
isdefined(@__MODULE__, :ref_itp) || (const ref_itp = water_itp.ref)

##################################################

# the time variable, shared by all objects and the derivative with respect to it
@variables t            
ddt = Differential(t)

# All the functions depending on wind and water vector fields.
# Note that `water_itp` and `wind_itp` must be loaded before using these.
v_x(x, y, t) = water_itp.u(x, y, t)
v_y(x, y, t) =  water_itp.v(x, y, t)
Dv_xDt(x, y, t) = MaterialDerivativeX(water_itp, x, y, t)
Dv_yDt(x, y, t) = MaterialDerivativeY(water_itp, x, y, t)
u_x(x, y, t, α) = (1 - α) * water_itp.u(x, y, t) + α * wind_itp.u(x, y, t)
u_y(x, y, t, α) = (1 - α) * water_itp.v(x, y, t) + α * wind_itp.v(x, y, t)
Du_xDt(x, y, t, α) = (1 - α) * MaterialDerivativeX(water_itp, x, y, t) + α * MaterialDerivativeX(wind_itp, x, y, t) 
Du_yDt(x, y, t, α) = (1 - α) * MaterialDerivativeY(water_itp, x, y, t) + α * MaterialDerivativeY(wind_itp, x, y, t) 
ω(x, y, t) = Vorticity(water_itp, x, y, t)

# Each function must be registered to allow symbolic computations
@register_symbolic v_x(x, y, t)
@register_symbolic v_y(x, y, t)
@register_symbolic Dv_xDt(x, y, t)
@register_symbolic Dv_yDt(x, y, t)
@register_symbolic u_x(x, y, t, α)
@register_symbolic u_y(x, y, t, α)
@register_symbolic Du_xDt(x, y, t, α)
@register_symbolic Du_yDt(x, y, t, α)
@register_symbolic ω(x, y, t)

# Symbolic functions computing the spring force
@register_symbolic spring_force_x(x1, x2, y1, y2, parameters)
@register_symbolic spring_force_y(x1, x2, y1, y2, parameters)

"""
    Clump(xy0, clump_parameters; name, forced)

Create an `ODESystem` representing a single clump which obey the BOM equations.

### Arguments

- `xy0` [km]: A vector with two elements representing the initial `[x, y]` coordinates of the clump.
- `clump_parameters`: A [`ClumpParameters`](@ref) object which contains the physics parameters of the clump.

### Optional Arguments

- `name`: A `Symbol` which is passed down to the name of the returned `ODESystem`.
- `forced`: A `Bool` which controls the extra force terms `Fx` and `Fy` in the equations. If `forced = false`, then these force terms are identically zero. If `forced = true` then the returned `ODESystem` does not contain equations for the `Fx` and `Fy` so that they can be provided manually.
"""
function Clump(
    xy0::Vector{<:Real},
    clump_parameters::ClumpParameters;
    name::Symbol,
    forced::Bool = false)

    @assert length(xy0) == 2 "Must provide `xy0` as a vector with two elements."

    ps = @parameters α = clump_parameters.α τ = clump_parameters.τ R = clump_parameters.R f = clump_parameters.f
    @variables x(t) = xy0[1] y(t) = xy0[2]
    @variables U(t)[1:2]
    @variables τU(t)[1:2]
    @variables τF(t)[1:2]

    eqs_U = [
        U[1] ~ u_x(x, y, t, α),
        U[2] ~ u_y(x, y, t, α)
    ]

    eqs_τU = [
        τU[1] ~ τ * (R*Dv_xDt(x, y, t) - R*(f + ω(x, y, t)/3)*v_y(x, y, t) - Du_xDt(x, y, t, α) + (f + R*ω(x, y, t)/3)*u_y(x, y, t, α)),
        τU[2] ~ τ * (R*Dv_yDt(x, y, t) + R*(f + ω(x, y, t)/3)*v_x(x, y, t) - Du_yDt(x, y, t, α) - (f + R*ω(x, y, t)/3)*u_x(x, y, t, α))
    ]

    eqs_xy = [
        ddt(x) ~ U[1] + τU[1] + τF[1],
        ddt(y) ~ U[2] + τU[2] + τF[2]
    ]

    eqs = [eqs_U ; eqs_τU ; eqs_xy]

    if !forced
        forces = [τF[1] ~ 0, τF[2] ~ 0]
        eqs = [eqs ; forces]
    end    


    return ODESystem(eqs, t, [x, y, U[1:2]..., τU[1:2]..., τF[1:2]...], ps; name)
end

"""
    Raft(xy0, clump_parameters, spring_parameters; name)
    

Create an `ODESystem` representing a network of `N` clumps connected by springs; the eBOM equations.

### Arguments

- `xy0` [km]: An `N x 2` matrix such that `xy0[i,:]` gives the `[x, y]` coordinates of the `i`th clump.
- `clump_parameters`: A vector of length `N` of [`ClumpParameters`](@ref) objects such that `clump_parameters[i]` are the parameters for the `i`th clump.
- `spring_parameters`: An `N x N` symmetric matrix of [`SpringParameters`](@ref) objects such that `spring_parameters[i, j]` gives the spring parameters for the spring connecting clump `i` with clump `j`. This matrix should be symmetric.

### Optional Arguments

- `name`: A `Symbol` which is passed down to the name of the returned `ODESystem`.
"""
function Raft(
    xy0::Matrix{<:Real},
    clump_parameters::Vector{<:ClumpParameters},
    spring_parameters::Matrix{<:SpringParameters};
    name::Symbol)

    @assert size(xy0, 2) == 2 "`xy0` must be an N x 2 matrix."
    @assert size(xy0, 1) == length(clump_parameters) == size(spring_parameters, 1) == size(spring_parameters, 2) "The dimensions of the input variables must be consistent."

    # a raft is a composition of several clumps
    N_clumps = size(xy0, 1)
    @named clump 1:N_clumps i -> Clump(xy0[i,:], clump_parameters[i], forced = true)

    # the force on each clump is the sum of the spring forces between it and every other clump
    forces_x = [clump[i].τF[1] ~ clump[i].τ*sum([spring_force_x(clump[i].x, clump[j].x, clump[i].y, clump[j].y, spring_parameters[i, j])] for j = 1:N_clumps if j != i)[1]
        for i = 1:N_clumps
    ]

    forces_y = [clump[i].τF[2] ~ clump[i].τ*sum([spring_force_y(clump[i].x, clump[j].x, clump[i].y, clump[j].y, spring_parameters[i, j])] for j = 1:N_clumps if j != i)[1]
        for i = 1:N_clumps
    ]

    # we keep track of the center of mass as an observable
    @variables COM(t)[1:2]
    eq_com = [COM[1] ~ sum([clump[i].x for i = 1:N_clumps])[1]/N_clumps, COM[2] ~ sum([clump[i].y for i = 1:N_clumps])[1]/N_clumps]

    eqs = [forces_x ; forces_y ; eq_com]

    return compose(ODESystem(eqs, t; name = name), clump[1:N_clumps]...)
end

"""
    RectangularRaft(x_range, y_range, clump_parameters, spring_parameters; network_type, name)

Construct a [`Raft`](@ref) in a rectangular arrangement such that each clump and spring have the same parameters.

### Arguments

- `x_range` [km]: A range which gives the x coordinates of the clumps in the raft. Should be increasing.
- `y_range` [km]: A range which gives the y coordinates of the clumps in the raft. Should be increasing.
- `clump_parameters`: The [`ClumpParameters`](@ref) shared by each clump.
- `spring_parameters`: The [`spring_parameters`](@ref) shared by each spring.

### Optional Arguments

- `network_type`: How the springs are conencted in the raft.
    - `"nearest"`: The default value. Each clump is connected to its perpendicular neighbors.
    - `"full"`: Each clump is connected to each other clump.
    - `"none"`: No clumps are connected.
- `name`: A `Symbol` which is passed down to the name of the returned `ODESystem`.
"""
function RectangularRaft(
    x_range::AbstractRange{<:Real}, 
    y_range::AbstractRange{<:Real},
    clump_parameters::ClumpParameters, 
    spring_parameters::SpringParameters; 
    network_type::String = "nearest",
    name::Symbol)

    @assert network_type in ["nearest", "full", "none"] "`network_type` not recognized."
    @assert step(x_range) > 0 "x range should be increasing."
    @assert step(y_range) > 0 "y range should be increasing."

    # a rectangular mesh
    network = reverse.(collect(Iterators.product(reverse(y_range), x_range))) # reverse so that the first row has the largest y
    n_col = length(x_range)
    n_row = length(y_range)
    N_clumps = length(network)

    xy0 = Matrix{Float64}(undef, N_clumps, 2)
    clump_parameters_raft = Vector{ClumpParameters}(undef, N_clumps)
    spring_parameters_raft = Matrix{SpringParameters}(undef, N_clumps, N_clumps)

    # these arrays are all constructed in dictionary order (across rows, then down columns)
    n(i, j) = (i - 1) * n_col + j

    for i = 1:n_row
        for j = 1:n_col
            # initial conditions
            xy0[n(i, j), :] .= network[i, j] 

            # all clumps have the same parameters
            clump_parameters_raft[n(i, j)] = clump_parameters

            # identify the appropriate connections
            if network_type == "nearest"
                connections = filter(idx -> (1 <= idx[1] <= n_row) && (1 <= idx[2] <= n_col), [(i-1, j), (i+1, j), (i, j-1), (i, j+1)])
            elseif network_type == "full"
                connections = [(a, b) for a = 1:n_row for b = 1:n_col]
            elseif network_type == "none"
                connections = []
            end
            
            connections = map(x->n(x...), connections)

            # all connected springs have the same parameters, `spring_parameters`
            # all unconnected springs have SpringParameters(k->0, 0) (effectively no spring)
            spring_parameters_raft[n(i, j), connections] .= spring_parameters
            spring_parameters_raft[n(i, j), setdiff(1:N_clumps, connections)] .= SpringParameters(k->0, 0)
        end
    end

    return Raft(xy0, clump_parameters_raft, spring_parameters_raft, name = name)
end


#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

# ref = EquirectangularReference(lon0 = -75.0, lat0 = 10.0);

# @load "water_itp.jld"
# @load "wind_itp.jld" 

# xy01 = [1.0, 2.0]
# xy0 = [1.0 2.0 ; 3.0 4.0]

# cp1 = ClumpParameters(ref);
# cps = [cp1, cp1]

# k_const(d) = 3
# sp = [SpringParameters(k_const, 0) SpringParameters(k_const, 4) ; SpringParameters(k_const, 6) SpringParameters(k_const, 0)];

# @named clump_no_force = Clump(xy01, cp1)
# @named clump_with_force = Clump(xy01, cp1, forced = true)
# @named raft = Raft(xy0, cps, sp)

# x_range = range(start = 5.0, length = 2, stop = 16.8)
# y_range = range(start = 3.0, length = 2, stop = 30.7)
# clump_parameters = ClumpParameters(ref)
# spring_parameters = SpringParameters(k -> 2.1, 1.5)

# @named rRaft = RectangularRaft(x_range, y_range, clump_parameters, spring_parameters)