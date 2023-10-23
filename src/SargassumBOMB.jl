include(joinpath(@__DIR__, "geography.jl"))
include(joinpath(@__DIR__, "physics.jl"))
include(joinpath(@__DIR__, "biology.jl"))
include(joinpath(@__DIR__, "trajectories.jl"))
include(joinpath(@__DIR__, "control.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/geo-methods.jl"))
include(joinpath(@__DIR__, "../../CustomMakie.jl/src/statistic-methods.jl"))

using SargassumFromAFAI

dists = SargassumDistribution(joinpath(@__DIR__, "..", "..", "SargassumFromAFAI.jl", "data", "dist-2018.nc"))

using Random: seed!

##############################################################################################################

