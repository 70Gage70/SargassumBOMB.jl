include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

using SargassumFromAFAI

isdefined(@__MODULE__, :DISTS_2018) || (const DISTS_2018 = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc")))
isdefined(@__MODULE__, :SFA_plot) || (const SFA_plot(time, week) = SargassumFromAFAI.plot(DISTS_2018[time], week, legend = false, resolution = (1920, 1080), limits = (-100, -40, 5, 35)))
isdefined(@__MODULE__, :SFA_plot!) || (const SFA_plot!(axis, time, week) = SargassumFromAFAI.plot!(axis, DISTS_2018[time], week))

using Random: seed!

##############################################################################################################

