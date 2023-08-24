using GeoDatasets
using JLD2

include(joinpath(@__DIR__, "..", "itp-core.jl"))

#############################################################################

@info "Constructing land interpolant."

# 0 is ocean, 1 is land and 2 is lake
lon, lat, land = GeoDatasets.landseamask(resolution = 'c', grid = 5)
land[land .== 2] .= 1 # lake is not ocean, so it's land

itp = GriddedField([:lon, :lat], Dict(:lon => lon, :lat => lat), Dict(:land => land), nothing, nothing, nothing, ref_itp)
itp = itp |> sph2xy |> x -> interpolate(x, interpolant_type = "nearest")

outfile = joinpath(@__DIR__, "land_itp.jld2")
rm(outfile, force = true)
jldsave(outfile, land_itp = itp)

@info "Land interpolant written to $(outfile)."





