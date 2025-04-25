# HDM_Connectivity_Analysis.R

## 1. Description
This script produces main figure outputs (e.g., richness, elevation, glaciation maps) and analyzes spatial correlations between residuals and predictors like cost distance across elevation. It uses the results of the moving window analysis to generate visualizations for the main manuscript figures.

## 2. Inputs (sorted by path)
**./File/**
- 1000m_UTM47_DEM.tif: Elevation raster (coarse resolution)
- 250m_UTM47.tif: Elevation raster (fine resolution)
- river.shp: River shapefile
- HDS_county.shp: County-level shapefile for the Hengduan Mountains
- end_div.tif: Endemic species richness raster
- glacial.tif: LGM glaciation raster
- AFEAD_HDM_faults.shp: Fault lines shapefile
- ele_summ.csv: Elevation-stratified summary statistics
- flux_mean_d10.tif: Flux raster
- Velo.tif: Climate velocity
- tem.tif: Temperature raster
- swb.tif: Soil water balance
- beta_cata.tif: Calculated endemic beta diversity categories
- tree1.rds: Dendrogram structure of beta diversity categories
- NMDS3_map.tif: NMDS3 analysis based on calculated species composition

**./Results/**
- res_rel_imp_average_50.tif: Moving window correlation result