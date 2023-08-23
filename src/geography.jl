using JLD2

include(joinpath(@__DIR__, "..", "interpolants", "interpolant-core.jl"))

########################################################################

# loading land interpolant and default EquirectangularReference
itp_path = joinpath(@__DIR__, "..", "interpolants", "land")
isdefined(@__MODULE__, :land_itp) || (const land_itp = load(joinpath(itp_path, "land_itp.jld2"), "land_itp"))
isdefined(@__MODULE__, :ref_itp) || (const ref_itp = land_itp.ref)

struct Land{I<:InterpolatedField, U<:Integer}
    land_itp::I
    deaths::Vector{U}
end

# condition 
function (land::Land)(u, t, integrator)
    land.deaths = [land.land_itp.fields[:land](clump_i(u, i)...) == 1.0 for i = 1:n_clumps(u)]
    return land.deaths > 0
end

# affect!
function (land::Land)(integrator)
    kill!(integrator, land.deaths)
end