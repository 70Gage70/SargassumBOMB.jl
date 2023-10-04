using NetCDF
using MAT
using Dates
using Unitful # ms = "m/s", val = 2.0 * uparse(ms)

####################################################################

raw_path = joinpath(@__DIR__, "..", "raw")
rick_raw = matopen(joinpath(raw_path, "merged_2018.mat"))

NetCDF_units = Dict(
    "m s**-1" => u"m/s",
    "m s-1" => u"m/s", 
    "degrees_C" => u"°C",
    "mmol m-3" => u"mmol/m^3"
)
const UNITS_OUT_SPEED = u"km/d"
const UNITS_OUT_TEMP = u"°C"
const UNITS_OUT_NUTR = u"mmol/m^3"

const t_start = DateTime(2018, 1, 1)
const t_end = DateTime(2018, 12, 31)

###########################
########################### RICK'S DATA
###########################

# indep vars: "Lat", "Lon", "T" (in that order!)
# dep vars: "Uek_a", "Uek_bar", "Ug_a", "Ug_bar", "Uslip_d", "Uslip_ud" and similarly for V
# longitude is in degrees E/W [-180, 180]
# latitude is in degrees N/S [-90 - 90]
# time is in Julian days
# V/U are in m/s

# to compute the actual water velocities:
# U = Uek_bar + Uek_a + Ug_bar + Ug_a + Uslip_d (ZONAL - x)
# V = Vek_bar + Vek_a + Vg_bar + Vg_a + Vslip_d (MERIDIONAL - y)

lon_rick = read(rick_raw, "Lon") |> vec .|> Float64
lon_idx = findall(x -> -101 <= x <= -49, lon_rick)
lon_rick = lon_rick[lon_idx]

lat_rick = read(rick_raw, "Lat") |> vec .|> Float64
lat_idx = findall(x -> 0 <= x <= 40, lat_rick)
lat_rick = lat_rick[lat_idx]

time_rick = read(rick_raw, "T") |> vec .|> julian2datetime .|> datetime2rata
t_start_idx = searchsortedfirst(rata2datetime.(time_rick), t_start)
t_end_idx = searchsortedfirst(rata2datetime.(time_rick), t_end)
time_rick = time_rick[t_start_idx:t_end_idx]

# u

Uek_bar = read(rick_raw, "Uek_bar")[lat_idx, lon_idx] |> x -> permutedims(x)
Uek_a = read(rick_raw, "Uek_a")[lat_idx, lon_idx, t_start_idx:t_end_idx] |> x -> permutedims(x, (2, 1, 3))
Uek =  zeros(eltype(Uek_a), size(Uek_a))
for t = 1:size(Uek_a, 3)
    Uek[:,:,t] = Uek_a[:,:,t] + Uek_bar
end

Ug_bar = read(rick_raw, "Ug_bar")[lat_idx, lon_idx] |> x -> permutedims(x)
Ug_a = read(rick_raw, "Ug_a")[lat_idx, lon_idx, t_start_idx:t_end_idx] |> x -> permutedims(x, (2, 1, 3))
Ug =  zeros(eltype(Ug_a), size(Ug_a))
for t = 1:size(Ug_a, 3)
    Ug[:,:,t] = Ug_a[:,:,t] + Ug_bar
end

Uslip_d = read(rick_raw, "Uslip_d")[lat_idx, lon_idx, t_start_idx:t_end_idx] |> x -> permutedims(x, (2, 1, 3))

u_rick = Uek + Ug  + Uslip_d

scale_factor = 1.0
units_factor = "m s-1" |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = NaN

u_rick[u_rick .== fill_value] .= 0.0
u_rick = u_rick .|> x -> scale_factor * units_factor * x 

# v

Vek_bar = read(rick_raw, "Vek_bar")[lat_idx, lon_idx] |> x -> permutedims(x)
Vek_a = read(rick_raw, "Vek_a")[lat_idx, lon_idx, t_start_idx:t_end_idx] |> x -> permutedims(x, (2, 1, 3))
Vek =  zeros(eltype(Vek_a), size(Vek_a))
for t = 1:size(Vek_a, 3)
    Vek[:,:,t] = Vek_a[:,:,t] + Vek_bar
end

Vg_bar = read(rick_raw, "Vg_bar")[lat_idx, lon_idx] |> x -> permutedims(x)
Vg_a = read(rick_raw, "Vg_a")[lat_idx, lon_idx, t_start_idx:t_end_idx] |> x -> permutedims(x, (2, 1, 3))
Vg =  zeros(eltype(Vg_a), size(Vg_a))
for t = 1:size(Vg_a, 3)
    Vg[:,:,t] = Vg_a[:,:,t] + Vg_bar
end

Vslip_d = read(rick_raw, "Vslip_d")[lat_idx, lon_idx, t_start_idx:t_end_idx] |> x -> permutedims(x, (2, 1, 3))

v_rick = Vek + Vg  + Vslip_d

scale_factor = 1.0
units_factor = "m s-1" |> x -> uconvert(UNITS_OUT_SPEED, 1.0 * NetCDF_units[x]).val
fill_value = NaN

v_rick[v_rick .== fill_value] .= 0.0
v_rick = v_rick .|> x -> scale_factor * units_factor * x 

####################################################################

rickout = joinpath(@__DIR__, "rick-2018.mat")

rm(rickout, force = true)

matwrite(rickout, Dict("lon" => lon_rick, "lat" => lat_rick, "t" => time_rick, "u" => u_rick, "v" => v_rick))