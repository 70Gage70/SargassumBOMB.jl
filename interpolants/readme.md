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

## Preprocessing 

The data requires preprossing to extract the relevant variables, ensure the units are correct and scale the time domain appropriately.  This is accomplished 
by `~/data/preprocessed/raw-preprocess.jl` The resulting data are stored in `.mat` files in `~/data/preprocessed`.

## Construction

There are three kinds of interpolants, `biology`, `land` and `ocean-atmos`, each stored in `~/interpolants/` in a folder of the same name. After the data 
are preprocessed, the actual interpolants may be constructed using the `itp-construct.jl` file in each folder. The interpolants are stored as [`.jld2`](https://github.com/JuliaIO/JLD2.jl) files.

- `biology`: Constructs `no3_itp.jld2` storing nitrogen and `temp_itp.jld2` storing temperature.
- `land`: Constructs `land_itp.jld2` storing a heatmap of land vs. ocean locations.
- `ocean-atmos`: Constructs `water_itp.jld2` storing ocean currents and `wind_itp.jld2` storing wind velocities.