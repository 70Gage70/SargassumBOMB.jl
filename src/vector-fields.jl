using DifferentialEquations
using ModelingToolkit
using Interpolations
using MAT

# xs = [x^2 for x = 1:0.2:5]
# ys = [y^1.2 for y = 1:0.1:5]
# A = [x^2 + y^3 for x in xs, y in ys]
# itp = linear_interpolation((xs, ys), A)

# itp(3.0, 3.0)
# gradient(itp, 3.0, 3.0)

"""
This is wind velocity data from 2021 in the north Atlantic. The dataset contains 5 keys:

"Lon" [399 x 1] [-99.75:0.25:-0.25 degrees]: The longitude (East/West) at which the measurement is taken.
"Lat" [360 x 1] [90.0:-0.25:0.25 degrees]: The latitude (North/South) at which the measurement is taken.
"t" [364 x 1] [738157.5:1.0:738520.5 Rata Die days}]: The time at which the measurement is taken. See `Dates.rata2datetime`.
"u" [360 x 399 x 364] [km/day]: The x component (East/West) of the wind velocity.
"v" [360 x 399 x 364] [km/day]: The y component (North/South) of the wind velocity.
"""

data_wind = matopen("../MATLAB/viento_2021.mat")
lon_wind, lat_wind, t_wind, u_wind, v_wind = read(data_wind, "Lon", "Lat", "t", "u", "v");
close(wind)

"""
This is water velocity data from 2021 in the north Atlantic. The dataset contains 5 keys:

"lon" [200 x 1] [-99.875:0.25:-50.125 degrees]: The longitude (East/West) at which the measurement is taken.
"lat" [120 x 1] [5.125:0.25:34.875 degrees]: The latitude (North/South) at which the measurement is taken.
"t" [365 x 1] [738157.0:1.0:738521.0 Rata Die days}]: The time at which the measurement is taken. See `Dates.rata2datetime`.
"u" [120 x 200 x 365] [km/day]: The x component (East/West) of the wind velocity.
"v" [120 x 200 x 365] [km/day]: The y component (North/South) of the wind velocity.
"""

data_wtr = matopen("../MATLAB/merged-2021-IAS.mat")
lon_wtr, lat_wtr, t_wtr, u_wtr, v_wtr = read(data_wtr, "lon", "lat", "t", "u", "v");
close(water)

"""
The dataset must now be cleaned. 
"""

# Reverse the order of the first dimension of u_wind and v_wind so that the latitudes are increasing.
u_wind = u_wind[end:-1:1, :, :]
v_wind = v_wind[end:-1:1, :, :]

# Set all NaNs (velocities measured over land) equal to zero.
u_wind[isnan.(u_wind)] .= 0.0
v_wind[isnan.(v_wind)] .= 0.0
u_wtr[isnan.(u_wtr)] .= 0.0
v_wtr[isnan.(v_wtr)] .= 0.0

# Rearrange u/v such that the order of entries is (lon, lat, t) == (x, y, t).
u_wind, v_wind, u_wtr, v_wtr  = [permutedims(x, [2, 1, 3]) for x in [u_wind, v_wind, u_wtr, v_wtr ]]

# We subtract the half-day from the wind velocity times and subtract the 365th day from the water velocity times and data,
# so that they are on the same grid.
# We also rescale the times to be in days since Jan. 1 2021.
t_wind = t_wind .- 0.5 .- 738157
t_wtr = t_wtr[1:end-1] .- 738157
u_wtr = u_wtr[:, :, 1:end-1]
v_wtr = v_wtr[:, :, 1:end-1]

# We convert lat, lon and t to iterators.
lon_wind, lat_wind, t_wind, lon_wtr, lat_wtr, t_wtr = [-99.75:0.25:-0.25, 0.25:0.25:90.0, 0.0:1.0:363.0, -99.875:0.25:-50.125, 5.125:0.25:34.875, 0.0:1.0:363.0]

"""
We construct the interpolators.
"""

bsp = BSpline(Cubic(Line(OnGrid()))) 

velocity_x_wind = scale(interpolate(u_wind, bsp), lon_wind, lat_wind, t_wind)
velocity_y_wind = scale(interpolate(v_wind, bsp), lon_wind, lat_wind, t_wind)
velocity_x_wtr = scale(interpolate(u_wtr, bsp), lon_wtr, lat_wtr, t_wtr)
velocity_y_wtr = scale(interpolate(v_wtr, bsp), lon_wtr, lat_wtr, t_wtr)
