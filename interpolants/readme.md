# Interpolants

## Introduction

The following data are required to run the *Sargassum* model:

- ocean currents (2d velocity field),
- ocean temperature (1d scalar field),
- 10 m wind velocities (2d velocity field).
- nutrients (scalar fields)

A dataset giving the location of landmass is also required, but this is already provided by [GeoDatasets.jl](https://github.com/JuliaGeo/GeoDatasets.jl).

Raw data should be uniformly gridded in space by longitude and latitude and time by day.

## Quickstart

### Downloading the raw data

Ensure that [`git-lfs`](https://git-lfs.com/) is installed on your system. Then, in the root directory of the package run
```
git lfs fetch
git lfs checkout
```

### Building the interpolants

In `~/interpolants`, run 
```julia
include("itp-construct-all.jl")
```

## Raw Data Sources

The data are sourced from the [Copernicus](https://www.copernicus.eu/en) Earth observation program, particularly the marine and climate subsets.
The data are stored in `~/data/raw` as NetCDF (`.nc`) files.

### Ocean Currents

Dataset: [GLOBAL_REANALYSIS_PHY_001_031](https://data.marine.copernicus.eu/product/GLOBAL_REANALYSIS_PHY_001_031/download)
Relevant variables: `uo_glor [m/s]`, `vo_glor [m/s]`

### Ocean Temperature

Dataset: [GLOBAL_REANALYSIS_PHY_001_031](https://data.marine.copernicus.eu/product/GLOBAL_REANALYSIS_PHY_001_031/download)
Relevant variables: `thetao_glor [°C]`

### Wind

Dataset: [ERA5](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview)
Relevant variables: `u10 [m/s]`, `v10 [m/s]`

### Nutrients

Dataset : [GLOBAL_MULTIYEAR_BGC_001_029](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_BGC_001_029/download)
Relevant variables: `no3 [mmol/m^3]`

### Waves and Stokes Drift

Dataset: [ERA5](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview)
Relevant variables: `Significant height of combined wind waves and swell [m]`, `U-component stokes drift  [m/s]`, `V-component stokes drift  [m/s]`

## Preprocessing 

The data requires preprossing to extract the relevant variables, ensure the units are correct and scale the time domain appropriately.  This is accomplished 
by `~/data/preprocessed/raw-preprocess.jl` The resulting data are stored in `.mat` files in `~/data/preprocessed`.

## Construction

There are four kinds of interpolants, `biology`, `land`, `ocean-atmos` and `waves`, each stored in `~/interpolants/` in a folder of the same name. After the data 
are preprocessed, the actual interpolants may be constructed using the `itp-construct.jl` file in each folder. The interpolants are stored as [`.jld2`](https://github.com/JuliaIO/JLD2.jl) files.

- `biology`: Constructs `no3_itp.jld2` storing nitrogen and `temp_itp.jld2` storing temperature.
- `land`: Constructs `land_itp.jld2` storing a heatmap of land vs. ocean locations.
- `ocean-atmos`: Constructs `water_itp.jld2` storing ocean currents and `wind_itp.jld2` storing wind velocities.
- `waves`: Constants `waves_itp.jld2` storing significant wave heights and `stokes_itp.jld2` storing Stokes drift velocities.

## Rick's Data

Data: data are stored in files with the naming convention “merged_YYYY.mat” where YYYY is the year, starting with 1993.  Variables contained within are:

(The following are identical for each year, and need only be loaded once)
- `Lon`, `Lat`: longitude and latitude grid points
- `Uek_bar`, `Vek_bar`: time mean zonal and meridional wind-driven velocities (m/s)
- `Ug_bar`, `Vg_bar`: time mean zonal and meridional geostrophic velocities (m/s)
 
(The following variables change in each yearly file) 
- `T`: time in Julian days.  2448989.00 = 1 January 1993, 0000 UTC.
- `Uek_a`, `Vek_a`: daily wind-driven current anomalies (m/s)
- `Ug_a`, `Vg_a`: daily geostrophic current anomalies (m/s)
- `Uslip_d`, `Vslip_d`: daily slip velocities for drogued drifter (m/s)
- `Uslip_ud`, `Vslip_ud`: daily slip velocities for undrogued drifters (m/s)
- `erH`: errors in sea level anomaly (m)
 
The total zonal speed of the water at 15m depth is `Uek_bar+Uek_a+Ug_bar+Ug_a`.  If the term `Uslip_d` is also added, then the zonal speed of a drogued drifter is reproduced.  If instead the term `Uslip_ud` is added, this matches as closely as possible the zonal speed of an undrogued drifter.

Here, `Uslip_d` and `Vslip_d` are used.