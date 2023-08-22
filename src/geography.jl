using JLD2

include(joinpath(@__DIR__, "..", "interpolants", "interpolant-constructors.jl"))

########################################################################

# loading land interpolant and default EquirectangularReference
itp_path = joinpath(@__DIR__, "..", "interpolants", "land")
isdefined(@__MODULE__, :land_itp) || (const land_itp = load(joinpath(itp_path, "land_itp.jld2"), "land_itp"))
isdefined(@__MODULE__, :ref_itp) || (const ref_itp = land_itp.ref)

