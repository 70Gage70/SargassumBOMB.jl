module Examples

using SargassumBOMB
using Dates

export generate_rp_example

function generate_rp_example()
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

end # module