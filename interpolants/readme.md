# Interpolants

## Introduction

The following data are required to run the *Sargassum* model:

- ocean currents (2d velocity field),
- ocean temperature (1d scalar field),
- 10 m wind velocities (2d velocity field).
- nutrients (scalar fields)

A dataset giving the location of landmass is also required, but this is already provided by [GeoDatasets.jl](https://github.com/JuliaGeo/GeoDatasets.jl).

Raw data should be uniformly gridded in space by longitude and latitude and time by day.

## Raw Data Sources

The data are sourced from the [Copernicus](https://www.copernicus.eu/en) Earth observation program, particularly the marine and climate subsets.

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
Relevant variables: `no3 [mmol/m3]`






