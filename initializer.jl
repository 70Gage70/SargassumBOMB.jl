# in root of BOM.jl
import Pkg
Pkg.activate(pwd())
cd("src")

include(joinpath(@__DIR__, "src/data-files.jl"))
include(joinpath(@__DIR__, "src/vf-refactor.jl"))

water_grid = VectorField2DGridSPH(water_file_default, lon_lat_time_order=[2, 1, 3]);
ref = EquirectangularReference(-75, 10);
water_itp_eqr = VectorField2DInterpolantEQR(water_grid, ref);