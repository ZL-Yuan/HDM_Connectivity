### =========================================================================
### Initialise system
### =========================================================================
# Start with house cleaning
rm(list = ls()); graphics.off()

library(terra)
library(tidyverse)
library(sp)
library(raster)

options(scipen=999)

# get working directory
wd <- getwd()

# set initial working directory
setwd(dirname(wd))

### =================================
### load in data for the project
### =================================

ele <- rast("./File/1000m_UTM47_DEM.tif") # elevation from DEM
Flux <- rast("./File/flux_mean_d10.tif") # landscape flux (connectivity)
velo <- rast("./File/Velo.tif") # Calculated climate change velocity (absolute value)
tem <- rast("./File/tem.tif") # Temperature in terms of mean temperature of the warmest quarter (MTWQ)
swb <- rast("./File/swb.tif") # soil water balance
end_div <- rast("./File/end_div.tif") # endemic richness (alpha richness)

### =============================================================== ###
### Moving window analysis for comparing endemic species richness 
### (alpha richness) to landscape level flux(connectivity) in 
### different starting position and window sizes 
### =============================================================== ###

# === 1. Prepare data frame from raster stack ===
flux_all <- data.frame(
  ele = values(ele),        # elevation from DEM
  end_div = values(end_div),# endemic richness (alpha richness)
  tem = values(tem),        # Temperature in terms of mean temperature of the warmest quarter (MTWQ)
  swb = values(swb),        # soil water balance
  Flux = values(Flux),      # landscape flux (connectivity)
  velo = values(velo) %>% abs() # Calculated climate change velocity (absolute value)
)

# Rename the columns for easier reference
names(flux_all) <- c("ele", "end_div", "tem", 'swb', 'Flux', 'velo')

# Add spatial coordinates (x, y) to the data frame
flux_all <- flux_all %>% 
  cbind(ele %>%
          xyFromCell(., 1:ncell(ele)) %>%
          as.data.frame()) %>%
  na.omit() %>%                    # remove rows with NA values
  filter(is.finite(Flux))         # keep only rows where Flux is finite

# === 2. Fit a Poisson GLM model to predict endemic richness ===
all_glm <- glm(
  end_div ~ tem + swb + log(velo) + I(tem^2) + I(swb^2) + I(log(velo)^2),
  family = "poisson",
  data = flux_all
) %>% step(trace = F) # Stepwise selection to simplify model

# Extract richness residuals and add them to the data frame
flux_all$resi <- residuals(all_glm)


# === 3. Initialize empty raster to hold moving window correlations ===
res_rel_imp <- raster(tem)
values(res_rel_imp) <- NA

# === 4. Define starting offsets to test different grid alignments ===
start_point <- list(
  c(0, 0),     # no offset
  c(0, 0.5),   # shift half a window up
  c(0.5, 0),   # shift half a window right
  c(0.5, 0.5)  # shift diagonally
)

# === 5. Loop through multiple spatial window sizes and offsets ===
for (window_size in c(50, 30, 20, 10)) {  # in kilometers
  for (i in 1:4) {  # loop through 4 offset start positions
    res_rel_imp <- raster(tem)
    values(res_rel_imp) <- NA
    
    # Slide window horizontally
    for (xmin in seq(-248500 + window_size*1000*start_point[[i]][1], 
                     (1123500 - window_size*1500), 
                     window_size*1000)) {
      cat(xmin, "\n")
      
      # Slide window vertically
      for (ymin in seq(2600000 + window_size*1000*start_point[[i]][2], 
                       (3803500 - window_size*1500), 
                       window_size*1000)) {
        cat(ymin, "\r")
        
        # Extract flux values within the current window
        flux_df <- flux_all[flux_all$x > (xmin + window_size*500) & flux_all$x < (xmin + window_size*1500),]
        flux_df <- flux_df[flux_df$y > (ymin + window_size*500) & flux_df$y < (ymin + window_size*1500),]
        
        flux_df <- na.omit(flux_df)
        flux_df <- flux_df[is.finite(flux_df$Flux),]
        
        # Limit analysis to mid-elevation data
        flux_df <- flux_df[flux_df$ele > 1000 & flux_df$ele < 5000, ] %>%
          na.omit()
        
        if (nrow(flux_df) < 50) {
          next # Skip if too few points in the window
        }
        
        # Define window extent and calculate correlation
        extent_window <- extent(c(xmin + window_size*500, xmin + window_size*1500, 
                                  ymin + window_size*500, ymin + window_size*1500))
        
        # Store Spearman correlation between residuals and flux
        res_rel_imp[extent_window] <- cor(flux_df$resi, flux_df$Flux, method = "spearman")
      }
    }
    
    # Save the result for this window size and offset
    res_rel_imp <- rast(res_rel_imp)
    writeRaster(res_rel_imp, file = paste0("./Results/res_rel_imp_", window_size, "_", i, ".tif"))
  }
}

# === 6. Average results across offsets for each window size ===
# 
for (window_size in c(50, 30, 20, 10)) {
  file_name <- paste0("./Results/res_rel_imp_", window_size, "_", 1:4, ".tif")
  res_rel_imp <- mean(rast(file_name), na.rm = TRUE) # stack & average
  writeRaster(res_rel_imp, file = paste0("./Results/res_rel_imp_average_", window_size, ".tif"))
}

### =============================================================== ###
### Moving window analysis for comparing endemic beta diversity 
### (beta richness) to landscape level flux(connectivity) in 
### different starting position and window sizes in focual area
### =============================================================== ###

# Load rasters
ele <- rast("./File/1000m_UTM47_DEM.tif") # elevation from DEM
Flux <- rast("./File/flux_mean_d10.tif") # landscape flux (connectivity)
velo <- rast("./File/Velo.tif") # Calculated climate change velocity (absolute value)
tem <- rast("./File/tem.tif") # Temperature in terms of mean temperature of the warmest quarter (MTWQ)
swb <- rast("./File/swb.tif") # soil water balance
beta <- rast("./File/beta_map.tif")             # Beta-diversity response variable

# Ensure all rasters are spatially aligned (same extent/resolution)
ele <- crop(ele, beta) %>% resample(beta)
Flux <- crop(Flux, beta) %>% resample(beta)
velo <- crop(velo, beta) %>% resample(beta)
tem <- crop(tem, beta) %>% resample(beta)
swb <- crop(swb, beta) %>% resample(beta)

# Create a unified data frame of raster values
flux_all <- data.frame(
  ele = values(ele),
  beta = values(beta),
  tem = values(tem),
  swb = values(swb),
  Flux = values(Flux),
  velo = values(velo) %>% abs()
)

# Name columns clearly
names(flux_all) <- c("ele", "beta", "tem", "swb", "Flux", "velo")

# Add X/Y coordinates and remove NA or non-finite Flux values
flux_all <- flux_all %>% 
  cbind(ele %>%
          xyFromCell(., 1:ncell(ele)) %>%
          as.data.frame()) %>%
  na.omit() %>% 
  filter(is.finite(Flux))

# Fit a GLM predicting beta-diversity using environmental variables and their squares
all_glm <- glm(
  beta ~ tem + swb + log(velo) + I(tem^2) + I(swb^2) + I(log(velo)^2),
  family = "poisson",
  data = flux_all
) %>% 
  step(trace = F)  # Stepwise selection to simplify model

# Save residuals from model to assess spatial correlation later
flux_all$resi <- residuals(all_glm)

# Prepare an empty raster to store correlation results
res_rel_imp <- raster(tem)
values(res_rel_imp) <- NA

# Define multiple starting points for offset grids (to reduce grid bias)
start_point <- list(
  c(0, 0),     # No offset
  c(0, 0.5),   # Offset vertically
  c(0.5, 0),   # Offset horizontally
  c(0.5, 0.5)  # Offset both directions
)

# === Moving window correlation analysis ===
# Loop through different window sizes (in km)
for (window_size in c(51, 30, 18, 12)) {
  for (i in 1:4) {  # loop through each offset setup
    res_rel_imp <- raster(tem)
    values(res_rel_imp) <- NA
    
    # Loop over X (longitude)
    for (xmin in seq(360500 + window_size * 1000 * start_point[[i]][1], 
                     (882500 - window_size * 1500), 
                     window_size * 1000)) {
      cat(xmin, "\n")
      
      # Loop over Y (latitude)
      for (ymin in seq(2849500 + window_size * 1000 * start_point[[i]][2], 
                       (3500500 - window_size * 1500), 
                       window_size * 1000)) {
        cat(ymin, "\r")
        
        # Subset flux_all to points inside the moving window
        flux_df <- flux_all[flux_all$x > (xmin + window_size * 500) & flux_all$x < (xmin + window_size * 1500),]
        flux_df <- flux_df[flux_df$y > (ymin + window_size * 500) & flux_df$y < (ymin + window_size * 1500),]
        
        flux_df <- na.omit(flux_df)
        flux_df <- flux_df[is.finite(flux_df$Flux),]
        
        # Focus only on mid-elevation band (removes high/low extremes)
        flux_df <- flux_df[flux_df$ele > 1000 & flux_df$ele < 5000, ] %>% na.omit()
        
        # Require at least 10 samples per window to compute correlation
        if (nrow(flux_df) < 10) {
          next
        }
        
        # Define the raster extent corresponding to this window
        extent_window <- extent(c(xmin + window_size * 500, xmin + window_size * 1500, 
                                  ymin + window_size * 500, ymin + window_size * 1500))
        
        # Calculate Spearman correlation between residuals and flux
        res_rel_imp[extent_window] <- cor(flux_df$resi, flux_df$Flux, method = "spearman")
      }
    }
    
    # Save the raster for this combination of window size and offset
    res_rel_imp <- rast(res_rel_imp)
    writeRaster(res_rel_imp, file = paste0("./Results/res_rel_imp_beta", window_size, "_", i, ".tif"), overwrite = TRUE)
  }
}

# === Average across the 4 offset rasters for each window size ===
for (window_size in c(51, 30, 18, 12)) {
  file_name <- paste0("./Results/res_rel_imp_beta", window_size, "_", 1:4, ".tif")
  
  # Load and average the 4 rasters (1 for each offset)
  res_rel_imp <- mean(rast(file_name), na.rm = TRUE)
  
  # Save the averaged map
  writeRaster(res_rel_imp, file = paste0("./Results/res_rel_imp_average_beta", window_size, ".tif"), overwrite = TRUE)
}