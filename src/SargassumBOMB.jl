include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

using SargassumFromAFAI

isdefined(@__MODULE__, :dists) || (const dists = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc")))
isdefined(@__MODULE__, :SFA_plot) || (const SFA_plot(time, week) = SargassumFromAFAI.plot(dists[time], week, legend = false, resolution = (1920, 1080), limits = (-100, -40, 5, 35)))

using Random: seed!

##############################################################################################################

