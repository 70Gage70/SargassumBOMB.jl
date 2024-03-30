using SargassumBOMB
using Test
using Dates

@testset "SargassumBOMB.jl" begin
	tspan = (DateTime(2018, 4, 13), DateTime(2018, 4, 15)) .|> datetime2time
	ics = InitialConditions(tspan, range(-55.0, -50.0, length = 5), range(5.0, 10.0, length = 5), to_xy = true)
    clumps = ClumpParameters()
    springs = BOMBSpring(1.0, ΔL(ics))
    connections = ConnectionsNearest(2)
    gd_model = BrooksModel(ics)
    land = Land()
    

    rp = RaftParameters(
        ics = ics,
        clumps = clumps,
        springs = springs,
        connections = connections,
        gd_model = gd_model,
        land = land
    )

    sol = simulate(rp, showprogress = false)
    @test true
end
