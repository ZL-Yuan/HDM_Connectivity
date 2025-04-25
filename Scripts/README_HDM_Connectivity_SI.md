# HDM_Connectivity_SI.R

## 1. Description
This script generates supplementary figures, including species occurrence maps, richness maps, predictor layers, and extended moving window visualizations for both alpha and beta diversity. It complements the main analysis by presenting additional insights.

## 2. Inputs (sorted by path)
**./File/**
- 1000m_UTM47_DEM.tif, 250m_UTM47.tif: Elevation data
- river.shp: River shapefile
- HDS_county.shp: County boundaries
- end_div.tif: Endemic richness
- glacial.tif: LGM glaciation raster
- AFEAD_HDM_faults.shp: Fault lines
- flux_mean_d10.tif: Flux raster
- Velo.tif: Climate change velocity raster
- tem.tif: MTWQ temperature
- swb.tif: Soil water balance

**./File_SI/**
- Allium_ovalifolium.csv: Species occurrence data
- allium ovalifolium.shp, allium ovalifolium.tif: Species range
- three_model_richness.tif: Overall richness map
- beta_map.tif: Beta diversity
- NMDS1_map.tif: NMDS1 analysis based on calculated species composition
- NMDS2_map.tif: NMDS2 analysis based on calculated species composition

**./Results/**
- res_rel_imp_average_10.tif, 20.tif, 30.tif, 50.tif: Moving window (alpha)
- res_rel_imp_average_beta12.tif, 18.tif, 30.tif, 51.tif: Moving window (beta)