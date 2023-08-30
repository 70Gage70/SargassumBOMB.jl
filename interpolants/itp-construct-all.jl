include(joinpath(@__DIR__, "..", "data", "preprocessed", "raw-preprocess.jl"))

include(joinpath(@__DIR__, "biology", "itp-construct.jl"))
include(joinpath(@__DIR__, "land", "itp-construct.jl"))
include(joinpath(@__DIR__, "ocean-atmos", "itp-construct.jl"))