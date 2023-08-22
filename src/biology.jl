using OrdinaryDiffEq
using JLD2

include(joinpath(@__DIR__, "parameters.jl"))
include(joinpath(@__DIR__, "..", "interpolants", "interpolant-constructors.jl"))

########################################################################

# loading interpolants
itp_path = joinpath(@__DIR__, "..", "interpolants", "ocean-atmos")
isdefined(@__MODULE__, :temp_itp) || (const temp_itp = load(joinpath(itp_path, "temp_itp.jld2"), "temp_itp"))
isdefined(@__MODULE__, :no3_itp) || (const no3_itp = load(joinpath(itp_path, "no3_itp.jld2"), "no3_itp"))