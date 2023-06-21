# set pwd here

include(joinpath(@__DIR__, "..", "src", "vector-fields", "vector-field-files.jl"))

construct_wind_itp_EQR()
construct_water_itp_EQR()
construct_itp_EQR()