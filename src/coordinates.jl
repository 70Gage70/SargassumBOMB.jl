"""
    sph2xy(lon, lat, lon0, lat0; R = 6371)

Compute planar coordinates (x, y) [km] from spherical coordinates (lon, lat) [deg] with references (lon0, lat0) [deg].

### Arguments

`lon`: Longitude in degrees (East/West).
`lat`: Latitude in degrees (North/South).
`lon0`: A reference longitude, usually taken near the middle of the domain.
`lat0`: A reference latitude, usually taken near the middle of the domain.
`R`: The radius of the Earth; x and y are returned in the units of R. Default 6371 km.
"""
function sph2xy(lon, lat, lon0, lat0; R = 6371)
    deg2rad = π/180

    x = R*(lon - lon0)*deg2rad*cos(lat0*deg2rad)
    y = R*(lat - lat0)*deg2rad

    return (x, y)
end

"""
    xy2sph(x, y, lon0, lat0; R = 6371)

Compute spherical coordinates (lon, lat) [deg] from planary coordinates (x, y) [km] from  with references (lon0, lat0) [deg].

### Arguments

`x`: The x Cartesian coordinate in km (East/West).
`y`: The y Cartesian coordinate in km (North/South).
`lon0`: A reference longitude, usually taken near the middle of the domain.
`lat0`: A reference latitude, usually taken near the middle of the domain.
`R`: The radius of the Earth; should be in the same units as x and y. Default 6371 km.
"""
function xy2sph(x, y, lon0, lat0; R = 6371)
    deg2rad = π/180
    rad2deg = 1/deg2rad

    lon = lon0 + rad2deg*x/(R*cos(lat0*deg2rad))
    lat = lat0 + rad2deg*y/R 

    return (lon, lat)
end