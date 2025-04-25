# HDM_moving_window.R

## 1. Description
This script performs moving window analyses to assess the spatial correlation between landscape connectivity (flux) and endemic species richness (both alpha and beta diversity) across varying window sizes and offsets in the Hengduan Mountains.

## 2. Inputs (sorted by path)
**./File/**
- 1000m_UTM47_DEM.tif: Elevation raster (DEM)
- flux_mean_d10.tif: Landscape connectivity (flux)
- Velo.tif: Climate change velocity
- tem.tif: Mean temperature of the warmest quarter (MTWQ)
- swb.tif: Soil water balance
- end_div.tif: Alpha (endemic) species richness
- beta_map.tif: Beta diversity map