# Project Overview

## General 
- This is the scripts used for the Connectivity project (*Yuan et al., 2025*)(DOI:   ). 
- The corresponding data is available at https://doi.org/10.5281/zenodo.15281421
- This file contains only the code needed for results and figure generation of the project, and it is to be used in combination with the code stated above.
- Please maintain the folder structure.


## Project Folder Structure
- `./Project/`
  - `Scripts/`: Contains the R scripts.
    - `HDM_moving_window.R`: Run this **first** to generate moving window correlation results.
    - `HDM_Connectivity_Analysis.R`: Run **second** to generate main figures and explore correlations.
    - `HDM_Connectivity_SI.R`: Run **last** for supplementary figures.
  - `File/`: Core data inputs for all scripts.
  - `File_SI/`: Additional data inputs used only in SI figures.
  - `Results/`: Stores outputs of moving window calculations.
  - `output/`: Final figure outputs for the main text.
  - `output_SI/`: Final figure outputs for the supplementary materials.


## Script Functions
- **HDM_moving_window.R**: Computes spatial correlations between endemic diversity residuals and landscape flux using moving window analysis.
- **HDM_Connectivity_Analysis.R**: Visualizes main analysis results, including figure panels and elevational relationships.
- **HDM_Connectivity_SI.R**: Generates extended analyses and supplementary visualizations.