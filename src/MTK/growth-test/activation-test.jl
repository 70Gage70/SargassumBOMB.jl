using ModelingToolkit, DifferentialEquations

struct LeafParameters{S<:Integer, T<:Real}
    x_speed::T
    y_speed::T
    activated::S
end

struct TreeParameters{T<:Real}
    stem_strength::T
end

@variables t
ddt = Differential(t)

function Leaf(xy0::Vector{<:Real}, lp::LeafParameters; name::Symbol)
    ps = @parameters act = lp.activated xsp = lp.x_speed ysp = lp.y_speed
    @variables x(t) = xy0[1] y(t) = xy0[2]
    @variables F(t)[1:2]

    eqs = [
        ddt(x) ~ act*(xsp + F[1]), 
        ddt(y) ~ act*(ysp + F[2])
    ]

    return ODESystem(eqs, t, [x, y, F[1:2]...], ps; name)
end

function Tree(xy0::Matrix{<:Real}, lp::Vector{<:LeafParameters}, tp::TreeParameters; name::Symbol)
    # a tree is a composition of leaves
    N_leaves = size(xy0, 1)
    @named leaf 1:N_leaves i -> Leaf(xy0[i,:], lp[i])
    ps = @parameters tps = tp.stem_strength

    # in practice, these would be complicated interaction forces
    forces_x = [
        leaf[i].F[1] ~ tps for i = 1:N_leaves
    ]

    forces_y = [
        leaf[i].F[2] ~ tps for i = 1:N_leaves
    ]

    eqs = [forces_x ; forces_y]

    return compose(ODESystem(eqs, t, [], ps; name = name), leaf[1:N_leaves]...)   
end

n_leaves = 3

xy0 = [rand() for i = 1:n_leaves, j = 1:2]
lp = [
    LeafParameters(1.1, 2.1, 1),
    LeafParameters(1.2, 2.2, 0),
    LeafParameters(1.3, 2.3, 0)
]

tp = TreeParameters(3.0)

@named tree = Tree(xy0, lp, tp)
tree = structural_simplify(tree)

t_range = (0.0, 10.0)

prob = ODEProblem(
    tree, 
    [],
    t_range
)

# prob.p = [3.0, 1.0, 1.1, 2.1, 0.0, 1.2, 2.2, 0.0, 1.3, 2.3]

function grow_leaf(t1::Real)
    condition(u, t, integrator) = t1 - t

    function affect!(integrator)
        u = integrator.u
        # resize!(integrator, length(u) + 2)
        u[3] = 0.0 # x
        u[4] = 0.0 # y

        integrator.p[5] = 1.0

        return nothing
    end

    return ContinuousCallback(condition, affect!)
end

@time sol = solve(prob, saveat = 1.0, callback = grow_leaf(3.0))