module Examples

using SargassumBOMB, SargassumFromAFAI
using Dates

export QuickRaftParameters

"""
    QuickRaftParameters()

Return a simple [`RaftParameters`](@ref) suitable for testing purposes.
"""
function QuickRaftParameters()
    LOADED_EXAMPLES = true
    try 
        itps_load(ITPS_DEFAULT_DIR)
    catch
        LOADED_EXAMPLES = false
    end

    tspan = (DateTime(2018, 4, 13), DateTime(2018, 4, 15)) .|> datetime2time
    ics = InitialConditions(tspan, range(-55.0, -50.0, length = 5), range(5.0, 10.0, length = 5), to_xy = true)
    n_clumps_max = size(ics.ics, 2)
    clumps = ClumpParameters()
    springs = BOMBSpring(1.0, ΔL(ics))
    connections = ConnectionsNearest(n_clumps_max, 2)
    gd_model = LOADED_EXAMPLES ? BrooksModel(n_clumps_max) : ImmortalModel(n_clumps_max)
    land = LOADED_EXAMPLES ? Land() : NoLand()

    return RaftParameters(
        ics = ics,
        clumps = clumps,
        springs = springs,
        connections = connections,
        gd_model = gd_model,
        land = land,
        n_clumps_max = n_clumps_max
    )
end


"""
    QuickRaftParameters(ymw_initial, ymw_final; kwargs...)

Generate a [`RaftParameters`](@ref) object to integrate from `(year, month, week)` inital to final. The raft 
is initialized at the Sargassum distribution at `ymw_initial`.

### Optional Arguments 

- `use_biology`: If `true`, include biological effects in the simulation. Default `false`.
- `n_levels`: This controls the number of initial clumps (higher is more clumps). \
See [`InitialConditions`](@ref) for more detail. Default `5`.
- `n_connections`: The number of inter-clump nearest neighbor connectionsto form. Default `2`. 
- `delta`: The bouyancy of the particle. Default: `1.14`.
- `tau`: The inertial response time, measured in days. Default `0.0103`.
- `sigma`: The Stokes drift coefficient. Default `1.20`.
- `A_spring`: The magnitude of the eBOMB spring stiffness. Default `15.1`.
- `lambda_spring`: A factor multiplying the springs' natural lengths. Default `2.97`.
- `mu_max`: Sargassum maximum growth rate, measured in inverse days. Default `0.00542`.
- `m`: Sargassum mortality rate, measured in inverse days. Default `0.00403`.
- `k_N`: Sargassum nutrient (N) uptake half saturation, measured in mmol N/m^3. Default `0.000129`.
- `S_min`: A clump dies when it's "amount" drops below this value. Default `-0.00481`.
- `S_max`: A clump dies when it's "amount" grows above this value. Default `0.001`.
- `seed`: A seed for reproducible randomness, passed to [`InitialConditions`](@ref). Default `1234`.
"""
function QuickRaftParameters(
	ymw_initial::NTuple{3, Integer},
	ymw_final::NTuple{3, Integer};
    use_biology::Bool = false,
	n_levels::Integer = 5, 
	n_connections::Integer = 2, 
	delta::Real = 1.14,
    tau::Real = 0.0103,
    sigma::Real = 1.20,
    A_spring::Real = 15.1,
    lambda_spring::Real = 2.97,
    mu_max::Real = 0.00542,
    m::Real = 0.00403,
    k_N::Real = 0.000129,
    S_min::Real = -0.00481,
    S_max::Real = 0.001,
	seed::Integer = 1234)
	
	tspan = (ymw2time(ymw_initial...), ymw2time(ymw_final...))
	dist_initial = SargassumFromAFAI.DIST_1718[ymw_initial[1:2]]
	ics = InitialConditions(tspan, dist_initial, [ymw_initial[3]], n_levels, seed = seed)
	n_clumps_max = size(ics.ics, 2)

	# Pick bio model
	if use_biology
		n_clumps_max_rp = 2*n_clumps_max
		gd_model = BrooksModel(
	        n_clumps_max_rp,
	        params = BrooksModelParameters(
	            μ_max = mu_max,
	            m = m,
	            k_N = k_N,
	            S_min = S_min,
	            S_max = S_max))
	else
		n_clumps_max_rp = n_clumps_max
		gd_model = ImmortalModel(n_clumps_max_rp)
	end

	clumps = begin c_p = ClumpParameters(δ = delta) ; ClumpParameters(c_p.α, tau, c_p.R, c_p.Ω, sigma) end

	connections = ConnectionsNearest(n_clumps_max_rp, n_connections)
	
	return RaftParameters(
	    ics = ics,
	    clumps = clumps,
	    springs = BOMBSpring(A_spring, lambda_spring*ΔL(ics)),
	    connections = connections,
	    gd_model = gd_model,
	    land = Land(),
	    n_clumps_max = n_clumps_max_rp,
	    fast_raft = false
	)
end

end # module