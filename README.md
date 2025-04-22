# README: Input Data for Code_v2_Figure1_Complete_Comments.R

This document describes the input datasets required for generating the figures using the script `code_v2_figure1_complete_comments.R`.

---

## Spatial Data Inputs

### 1. County Boundary Shapefile
- **File:** `./File/HDS_county.shp`
- **Description:** County-level shapefile for the Hengduan Mountains region. Used to define the study area.

### 2. Digital Elevation Models (DEM)
- **File:** `./File/250m_UTM47.tif`
  - High-resolution DEM (250m) used for hillshade and fine-scale terrain analysis.
- **File:** `./File/1000m_UTM47_DEM.tif`
  - Coarser DEM (1000m) used in general spatial analysis.

### 3. River Network
- **File:** `./File/river.shp`
- **Description:** Vector shapefile of river networks, projected and masked to the study area.

### 4. Endemic Richness Raster
- **File:** `./File/end_div.tif`
- **Description:** Raster showing endemic species richness, used in multiple figure panels.

### 5. Glacial Extent (Last Glacial Maximum)
- **File:** `./File/glacial.tif`
- **Description:** Raster showing areas covered by glaciers during the LGM. Used for visual overlays and comparison.

### 6. Tectonic Fault Lines
- **File:** `./File/AFEAD_HDM_faults.shp`
- **Description:** Vector file with geological fault lines for tectonic context.

---

## Tabular Data Inputs

### 7. Elevational Summary Table
- **File:** `./File/ele_summ.csv`
- **Description:** Elevation-band summary data including area, richness, and predictor values. Used for model fitting and plotting (Figures 1B, 2B, 2C).

---

## Analytical Results and Intermediate Inputs

### 8. Moving Window Correlation Results
- **File:** `./Results/res_rel_imp_average_50.tif`
  - Richness vs flux correlation (used in Figure 2A).
- **File:** `./Results/res_rel_imp_average_beta30.tif`
  - Beta diversity vs flux correlation (used in Figure 3C).

---

## Beta Diversity and Composition Analysis

### 9. Categorical Beta Diversity Map
- **File:** `./File/beta_cata.tif`
- **Description:** Raster of classified beta-diversity categories (12 types), used in Figures 3A and 3B.

### 10. Composition Tree (Dendrogram)
- **File:** `./File/tree1.rds`
- **Description:** RDS object of hierarchical clustering tree used for plotting group composition (Figure 3AB).

### 11. NMDS Ordination Map
- **File:** `./File/NMDS3_map.tif`
- **Description:** NMDS Axis 3 raster map representing compositional gradients (Figure 3D).

---

## Notes
- All files must be placed in the `./File/` or `./Results/` directories as specified.
- Output figures are saved to the `output/` directory in both PNG and PDF formats.

---

If you need help generating or locating these files, or want a data download script, feel free to ask!
