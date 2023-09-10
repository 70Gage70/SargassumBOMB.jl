using NetCDF

afai = "afai.nc"
ncread(afai, "longitude")

###########################
########################### WIND
###########################

# variables: longitude, latitude, time, AFAI(longitude, latitude, time)
# longitude is in degrees E/W [-180, 180]
# latitude is in degrees N/S [-90, 90]
# time is in seconds since 1970-01-01T00:00:00
# AFAI is an index