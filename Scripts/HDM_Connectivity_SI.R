### =========================================================================
### Initialise system
### =========================================================================

# Start with house cleaning
rm(list = ls()); graphics.off()  # Clear all objects and close plots

# Load required libraries
library(terra)         # spatial raster/vector handling
library(RColorBrewer)  # color palettes
library(magick)        # for image processing (e.g., combining plots)
library(tidyverse)     # tidy data manipulation
library(viridis)       # perceptually uniform color maps
library(sp)            # spatial objects (used by some older code)
library(raster)        # legacy raster package (some tools not yet in terra)
library(ggplot2)       # grammar of graphics
library(fields)        # image.plot() for legends
library(magrittr)

options(scipen = 999)  # avoid scientific notation in axis labels

# Get current working directory
wd <- getwd()

# Move one directory up
setwd(dirname(wd))

### =================================
### load in data for the project
### =================================

# Hengduan Mountains county-level shapefile
Hds_county <- vect('./File/HDS_county.shp') 
Hds_county_outline <- aggregate(Hds_county)  # dissolve boundaries to get single polygon

# Load DEM at 250m resolution
ele_raw <- rast("./File/250m_UTM47.tif")
ele_raw <- crop(ele_raw, Hds_county_outline)  # crop to study area
crs(ele_raw) <- CRS("+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs +type=crs")  # set UTM projection

# Load DEM at 1000m resolution for coarser analyses
ele <- rast("./File/1000m_UTM47_DEM.tif")

# Load river shapefile and prepare
river <- vect("./File/river.shp")
river <- terra::project(river, "+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs +type=crs")
river <- crop(river, Hds_county_outline)
river <- mask(river, Hds_county_outline)

# Load endemic richness raster and crop/mask to study area
end_div <- rast("./File/end_div.tif")
end_div <- crop(end_div, Hds_county_outline)
end_div <- mask(end_div, Hds_county_outline)

# Load glacial extent from LGM (Last Glacial Maximum)
glacial <- rast("./File/glacial.tif")
glacial[glacial == 0] <- NA  # turn 0 into NA to exclude unglaciated areas
glacial <- crop(glacial, Hds_county_outline)
glacial <- mask(glacial, Hds_county_outline)

# Load faults shapefile (tectonic features)
faults <- vect("./File/AFEAD_HDM_faults.shp")

### =========================================================================
### Prepare decorations
### =========================================================================

# Calculate slope and aspect for hillshading
slp <- terra::terrain(ele_raw, v = "slope", unit = "radians")
asp <- terrain(ele_raw, v = "aspect", unit = "radians")
hillsh <- shade(slp, asp, direction = 315)  # shaded relief map

# Semi-transparent black overlay for hillshade (used to darken hillshade)
hillcol <- colorRampPalette(c("#00000000", "#000000DA"), alpha = TRUE)

# Define color palette for endemic richness (multi-color gradient)
col_vec <- colorRampPalette(c("#9ECAE1", "#0099CC", "#58BC5D", "#EEF559", "#FF9933", "red"))(100)

# Custom ETH Zurich blue colors (used in line and polygon plots)
ETH_blue_old <- "#1F407A"
ETH_blue_new <- "#215CAF"

### =========================================================================
### Plot Allium ovalifolium occurrence and range Supplementary Figure 2
### 1) Occurrences only
### 2) Occurrences + range ploygon, highlight outliers
### 3) Range only
### =========================================================================

# Load and prepare occurrence data
occ <- read.csv2('./File_SI/Allium_ovalifolium.csv', sep = ',')
occ <- occ[, c(3, 2)]
colnames(occ) <- c('x', 'y')
occ <- occ[!apply(is.na(occ), 1, all), ]
occ$x <- as.numeric(occ$x)
occ$y <- as.numeric(occ$y)

# Convert to spatial points and project
spatial_data <- vect(occ, geom = c("x", "y"), crs = "EPSG:4326")
spatial_data_projected <- project(spatial_data, 'EPSG:32647')

# Load polygon range and raster
poly.spp <- vect('./File_SI/allium ovalifolium.shp')
poly.rast <- rast('./File_SI/allium ovalifolium.tif')

# Identify occurrence points outside the polygon using vector-based intersection
in_poly <- terra::relate(spatial_data_projected, poly.spp, "intersects")

# Points not intersecting the polygon
outliers <- spatial_data_projected[!in_poly, ]

# Function to set up plot layout
init_plot <- function(file_name) {
  pdf(file_name, , width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)
  par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), ps = 8, cex = 1)
  plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8)
  plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)
  plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 1.5)
}

# ====================
# Plot 1: Occurrences only
# ====================
init_plot("output_SI/sf2_plot1_occurrences_only.pdf")
plot(spatial_data_projected, add = TRUE, col = "black", pch = 16, cex = 1)
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = "a", font = 1, cex = 2)
legend(-150000, 2900000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 0.9)
dev.off()

# ====================
# Plot 2: Occurrences + polygon + red circles around outliers
# ====================
init_plot("output_SI/sf2_plot2_occurrences_and_range.pdf")
plot(poly.spp, add = TRUE, border = "#355E3B", col = "#355E3B80")
plot(spatial_data_projected, add = TRUE, col = "black", pch = 16, cex = 1)
points(crds(outliers), pch = 21, bg = NA, col = "red", lwd = 2, cex = 2.5)
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = "b", font = 1, cex = 2)
legend(-150000, 2900000, legend = c("Rivers", "Outliers"),
       col = c(ETH_blue_old, "red"), pch = c(NA, 21), lty = c(1, NA), lwd = c(2, 1.5),
       pt.bg = c(NA, NA), box.lty = 0, cex = 0.9)
dev.off()

# ====================
# Plot 3: Polygon range only
# ====================
init_plot("output_SI/sf2_plot3_range_only.pdf")
plot(poly.rast, add = TRUE, col = "#355E3B", alpha = 0.8, legend = FALSE)
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = "c", font = 1, cex = 2)
legend(-150000, 2900000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 0.9)
dev.off()

# ====================
# combine the plots
# ====================

# Define desired physical size (in cm) and resolution (in DPI)
width_cm <- 7.95
height_cm <- 6.9
dpi <- 300

# Convert dimensions to pixels
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

# Step 1: Read in the PDF plots (first page only)
fig1 <- image_read_pdf("output_SI/sf2_plot1_occurrences_only.pdf", density = dpi)
fig2 <- image_read_pdf("output_SI/sf2_plot2_occurrences_and_range.pdf", density = dpi)
fig3 <- image_read_pdf("output_SI/sf2_plot3_range_only.pdf", density = dpi)

# Step 2: Resize each image to the exact width and height
fig1 <- image_resize(fig1, paste0(target_width, "x", target_height, "!"))
fig2 <- image_resize(fig2, paste0(target_width, "x", target_height, "!"))
fig3 <- image_resize(fig3, paste0(target_width, "x", target_height, "!"))

# Step 3: Stack them vertically
combined_fig <- image_append(c(fig1, fig2, fig3), stack = TRUE)

# Step 4: Save the combined figure
image_write(combined_fig, path = "output_SI/sf2_combined_species_plots.png", format = "png")
image_write(combined_fig, path = "output_SI/sf2_combined_species_plots.pdf", format = "pdf")


### =========================================================================
### Plot overall richness Supplementary Figure 3
### =========================================================================

# Load overall richness raster
overall_richness <- rast('./File_SI/three_model_richness.tif')

# Define legend tick positions and labels
rn <- c(0, 2742)
lab_vals <- c(0, 500, 1000, 1500, 2000, 2500, 3000)
lab_vals <- lab_vals[lab_vals >= rn[1] & lab_vals <= rn[2]]
lab_names <- as.character(lab_vals)

# Start PNG device
# png("output/overall_richness.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)
pdf("output_SI/sf8_overall_richness.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)  # for publication-ready PDF

# Set plot parameters
par(mfrow=c(1,1),oma=c(0,0,0,0),ps=8,cex=1)

# Plot the richness raster without axes or legend
plot(overall_richness, col = col_vec, legend = FALSE, axes = FALSE, mar = c(0.1, 0.1, 0.1, 4.2), border = FALSE)

# Overlay hillshade for terrain relief
plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)

# Overlay county outlines
plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)

# Overlay river layer
plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 1.5)

# Define color breaks aligned with col_vec
n_col <- length(col_vec)
breaks <- seq(rn[1], rn[2], length.out = n_col + 1)

# Add vertical legend aligned with breaks and ticks
fields::image.plot(legend.only = TRUE,
                   zlim = rn,
                   col = col_vec,
                   breaks = breaks,
                   smallplot = c(0.87, 0.9, 0.2, 0.9),  # position of legend in figure space
                   legend.lab = "Overall richness",
                   axis.args = list(at = lab_vals, labels = lab_names, cex.axis = 0.8),
                   legend.line = -2)


# Add river legend
legend(-150000, 2900000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 0.9)

# Close graphics device
dev.off()

# Clean up memory
gc()

# Define desired physical size (in cm) and resolution (in DPI)
width_cm <- 7.95
height_cm <- 6.9
dpi <- 300

# Convert dimensions to pixels
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

# Step 1: Read in the first page of the PDF at high resolution
fig <- image_read_pdf("output_SI/sf3_overall_richness.pdf", density = dpi)

# Step 2: Resize image to exact pixel dimensions
fig <- image_resize(fig, paste0(target_width, "x", target_height, "!"))  # "!" forces exact size

# Step 3: Export as PNG
image_write(fig, path = "output_SI/sf3_overall_richness.png", format = "png")

### =========================================================================
###  Supplementary Figure 7 beta diversity at focal study region
### =========================================================================

# Read in beta diversity map
r <- raster('File_SI/beta_map.tif')

beta_scale <- r %>% scale()

# Standardize (scale) the beta diversity raster
beta_scaled <- (beta_scale - cellStats(beta_scale, stat='min', na.rm=TRUE)) / 
  (cellStats(beta_scale, stat='max', na.rm=TRUE) - cellStats(beta_scale, stat='min', na.rm=TRUE))

beta_scaled <- rast(beta_scaled)

# Define color scale limits and corresponding labeled ticks
rn <- c(0, 1)
lab_vals <- c(0,0.25,0.5,0.75,1)
lab_names <- as.character(lab_vals)

# Define color palette for visualization
beta_colors <- unique(colorRampPalette(c("steelblue1", "khaki1", "orangered"))(85))

# Compute breaks for the color legend to align with the color vector
n_col <- length(beta_colors)
breaks <- seq(rn[1], rn[2], length.out = n_col + 1)

# Start PNG device for output
#png("output/supp4.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)

# Optional PDF version
pdf("output_SI/sf7_beta_map.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

# Set up base plotting parameters
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Plot the scaled beta diversity raster without legend or axes
plot(beta_scaled, 
     range = rn, 
     col = beta_colors,
     plg = list(title = "Beta Diversity"),
     mar = c(0.5, 0.5, 0.1, 1), axes = FALSE, border = FALSE, alpha = 1, legend = FALSE)

# Overlay river lines
plot(river, col = ETH_blue_old, add = T, border = FALSE, lwd = 1)

# Add hillshade for terrain context
plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)

# Overlay county boundaries
plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)

# Add vertical legend using fields::image.plot()
par(xpd = NA)
# Add river legend with transparent background
legend(280000, 2980000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 1, bg = "transparent")
fields::image.plot(legend.only = TRUE,
                   zlim = rn,
                   col = beta_colors,
                   breaks = breaks,
                   smallplot = c(0.87, 0.9, 0.2, 0.9),  # position of the legend
                   legend.lab = "Beta Diversity",
                   legend.cex = 0.8,
                   axis.args = list(
                     at = lab_vals,
                     labels = lab_names,
                     cex.axis = 0.8
                   ),
                   legend.line = -1.8)



par(xpd = FALSE)

# Close graphics device
dev.off()
gc()

### =========================================================================
### Plot Supplementary Figure 8 - Predictor rasters (Flux, Velocity, MTWQ, SWB)
### =========================================================================

library(magrittr)
library(Cairo)

# Load rasters
Flux <- rast("./File/flux_mean_d10.tif") %>% mask(Hds_county_outline)
velo <- rast("./File/Velo.tif") %>% mask(Hds_county_outline)
tem <- rast("./File/tem.tif") %>% mask(Hds_county_outline)
swb <- rast("./File/swb.tif") %>% mask(Hds_county_outline)

# Color palette (same across maps)
pred_colors <- unique(colorRampPalette(c("steelblue1", "khaki1", "orangered"))(85))

# List of rasters and titles
rasters <- list(Flux, velo, tem, swb)
x <- expression("Climate Change Velocity (" * Delta * "m/yr)")
titles <- c("Flux", x, "Temperature (MTWQ, °C)", "Site Water Balance (mm/yr)")
titles_small <- c("Flux", "Climate Change Velocity", "Temperature", "Site Water Balance")
pic_order <- c('a','b','c','d')

# Loop through each raster layer
for (i in seq_along(rasters)) {
  # Open PDF output
  CairoPDF(str_c("output_SI/sf8_predictor_",titles_small[i],".pdf"), width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)
  par(mfrow=c(1,1),oma=c(0,0,0,0),ps=8,cex=1)
  r <- rasters[[i]]
  title_i <- titles[i]
  # Calculate data range and breaks automatically
  rn <- range(values(r), na.rm = TRUE)
  n_col <- length(pred_colors)
  breaks <- seq(rn[1], rn[2], length.out = n_col + 1)
  
  # Step 1: Create 5 evenly spaced values inside rn (excluding ends)
  inner_seq <- seq(rn[1], rn[2], length.out = 7)[2:6]
  
  # Step 2: Round up to nearest 5
  if(rn[2]>=5){
  rounded_seq <- ceiling(inner_seq/ 5) * 5
  clipped_inner <- pmin(pmax(rounded_seq, rn[1]), rn[2])
  sort_seq <- sort(unique(c(rn[1], clipped_inner, rn[2])))
  full_seq <- round(sort_seq, 0)
  }else{
    rounded_seq <- inner_seq
    clipped_inner <- pmin(pmax(rounded_seq, rn[1]), rn[2])
    sort_seq <- sort(unique(c(rn[1], clipped_inner, rn[2])))
    full_seq <- round(sort_seq, 1)
  }
  
  # Plot raster
  plot(r, col = pred_colors, legend = FALSE, axes = FALSE, mar = c(0.1, 0.1, 0.1, 4.2))
  
  # Overlay county outlines
  plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)
  
  # Overlay hillshade
  plot(hillsh, bty = "n", legend = FALSE, axes = FALSE,col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)
  
  par(xpd = NA)
  text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
       y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
       labels = pic_order[i], font = 1, cex = 2)
  par(xpd = FALSE)
  
  #plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)
  plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 1.5)
  
  # Add river legend
  legend(-150000, 2900000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 0.9)
  
  # Add vertical color legend for each map
  fields::image.plot(legend.only = TRUE,
                     zlim = rn,
                     col = pred_colors,
                     breaks = breaks,
                     smallplot = c(0.88, 0.9, 0.2, 0.9),  # right margin
                     legend.lab = title_i,
                     legend.cex = 0.8,
                     axis.args = list(at = sort_seq, labels = as.character(full_seq), cex.axis = 0.8),
                     legend.line = -1.8)
  
  # Close device
  dev.off()
  gc()
  
  print(paste0(title_i, ' DONE'))
}

# Define desired physical size (in cm) and resolution (DPI)
width_cm <- 7.95
height_cm <- 6.9
dpi <- 300

# Convert dimensions to pixels
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

# Read PDFs and convert to images
fig1 <- image_read_pdf("output_SI/sf8_predictor_Flux.pdf", density = dpi)
fig2 <- image_read_pdf("output_SI/sf8_predictor_Climate Change Velocity.pdf", density = dpi)
fig3 <- image_read_pdf("output_SI/sf8_predictor_Temperature (MTWQ).pdf", density = dpi)
fig4 <- image_read_pdf("output_SI/sf8_predictor_Site Water Balance.pdf", density = dpi)

# Resize all to same dimensions
fig1 <- image_resize(fig1, paste0(target_width, "x", target_height, "!"))
fig2 <- image_resize(fig2, paste0(target_width, "x", target_height, "!"))
fig3 <- image_resize(fig3, paste0(target_width, "x", target_height, "!"))
fig4 <- image_resize(fig4, paste0(target_width, "x", target_height, "!"))

# Stack images into a 2x2 grid
row1 <- image_append(c(fig1, fig2))
row2 <- image_append(c(fig3, fig4))
final_plot <- image_append(c(row1, row2), stack = TRUE)

# Save output
image_write(final_plot, path = "output_SI/sf8_predictor_combined.png", format = "png")
image_write(final_plot, path = "output_SI/sf8_predictor_combined.pdf", format = "pdf")

### =========================================================================
### Plot Supplementary Figure 9 - correlation between flux and endemic plant 
### species richness (α-diversity) residuals in different window size
### =========================================================================

## ===== labels ====
rn <- c(-1, 1)  # range of correlation values

# Legend tick marks and labels
lab <- c('-1' = -1,
         '-0.5' = -0.5,
         '0' = 0,
         '0.5' = 0.5,
         '1' = 1)
## =================

# Define diverging color palette (red to blue)
contrast.color <- colorRampPalette(rev(c(
  '#67001f',
  '#b2182b',
  '#d6604d',
  '#f4a582',
  '#fddbc7',
  '#F2F2F2',
  '#d1e5f0',
  '#92c5de',
  '#4393c3',
  '#2166ac',
  '#053061')))

# read in the moving window results file
res_rel_imp_50 <- rast('./Results/res_rel_imp_average_50.tif')
res_rel_imp_30 <- rast('./Results/res_rel_imp_average_30.tif')
res_rel_imp_20 <- rast('./Results/res_rel_imp_average_20.tif')
res_rel_imp_10 <- rast('./Results/res_rel_imp_average_10.tif')

rasters <- list(res_rel_imp_10,res_rel_imp_20,res_rel_imp_30,res_rel_imp_50)
titles <- c('10km','20km','30km', '50km')

#png("output/richvsflux_f2a.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)

for (i in 1:4){

pdf(str_c("output_SI/sf9_moving_window_",titles[i],".pdf"), width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

par(mar = c(0.1, 0.1, 0.1, 4.2), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

r <- rasters[[i]]
# Mask correlation raster with study boundary
res_masked <- mask(r, Hds_county_outline)

# Plot raster image
plot(res_masked, 
     range = c(-1, 1), 
     col = contrast.color(16),
     plg = list(title = ""),
     mar=c(0.5,0.1,0.1,4.2), axes=F, border=T , alpha = 0.7, legend = F)

# Add hillshade on top for terrain effect
plot(hillsh,bty="n",legend=F,axes=F,col=rev(hillcol(50)),maxcell=1e8,add=T)

# Add county outline
plot(Hds_county_outline, add = TRUE, lwd = 1, col = NA)

# Add river layer
plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 1.5)

# Add color legend using fields::image.plot

# Define break points for color edges (ensure exact alignment with color bins)
breaks <- seq(rn[1], rn[2], length.out = 17)  # 16 color bins = 17 breaks

fields::image.plot(legend.only = TRUE,
                   zlim = rn,
                   col = contrast.color(16),
                   breaks = breaks,  # enforce correct bin edges
                   smallplot = c(0.87, 0.9, 0.2, 0.9),  # legend position
                   legend.lab = "Spearman's Rank Correlation",
                   legend.cex = 0.8,
                   axis.args = list(
                     at = lab,                     # tick positions
                     labels = names(lab),          # display as "-1", "-0.5", ..., "1"
                     cex.axis = 0.8
                   ),
                   legend.line = -1.8)

# Add panel label 'a'
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = pic_order[i], font = 1, cex = 2)

# Add river legend
legend(-150000, 2800000,
       legend = "Rivers",
       col = ETH_blue_old,
       lty = 1, lwd = 2,
       box.lty = 0, cex = 0.9)

dev.off()
gc()
}

# Combining the figures

# Set image size and DPI
width_cm <- 7.95
height_cm <- 6.9
dpi <- 300
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

# Read each individual plot
fig1 <- image_read_pdf("output_SI/sf9_moving_window_10km.pdf", density = dpi)
fig2 <- image_read_pdf("output_SI/sf9_moving_window_20km.pdf", density = dpi)
fig3 <- image_read_pdf("output_SI/sf9_moving_window_30km.pdf", density = dpi)
fig4 <- image_read_pdf("output_SI/sf9_moving_window_50km.pdf", density = dpi)

# Resize to uniform dimensions
fig1 <- image_resize(fig1, paste0(target_width, "x", target_height, "!"))
fig2 <- image_resize(fig2, paste0(target_width, "x", target_height, "!"))
fig3 <- image_resize(fig3, paste0(target_width, "x", target_height, "!"))
fig4 <- image_resize(fig4, paste0(target_width, "x", target_height, "!"))

# Combine into 2x2 grid
row1 <- image_append(c(fig1, fig2))
row2 <- image_append(c(fig3, fig4))
combined <- image_append(c(row1, row2), stack = TRUE)

# Save combined image
image_write(combined, path = "output_SI/sf9_moving_window_combined.png", format = "png")
image_write(combined, path = "output_SI/sf9_moving_window_combined.pdf", format = "pdf")

### =========================================================================
### Plot Supplementary Figure 10 - correlation between flux and endemic plant 
### species richness (b-diversity) residuals in different window size
### =========================================================================

# read in moving window result of beta diversity

res_rel_imp_beta_51 <- rast('./Results/res_rel_imp_average_beta51.tif')
res_rel_imp_beta_30 <- rast('./Results/res_rel_imp_average_beta30.tif')
res_rel_imp_beta_18 <- rast('./Results/res_rel_imp_average_beta18.tif')
res_rel_imp_beta_12 <- rast('./Results/res_rel_imp_average_beta12.tif')

rasters <- list(res_rel_imp_beta_12,res_rel_imp_beta_18,res_rel_imp_beta_30,res_rel_imp_beta_51)
titles <- c('12km','18km','30km', '51km')

# Define color scale range
rn = c(-1, 1)

# Define breaks and labels for the color legend
lab = c('-1' = -1,
        "-0.5" = -0.5,
        "0" = 0,
        "0.5" = 0.5,
        "1" = 1)


#png("output/betavsflux_f3c.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)

for (i in 1:4){

pdf(str_c("output_SI/sf10_beta_window_",titles[i],".pdf"), width = 7.95/2.54, height = 6.9/2.54,pointsize=7.5, useDingbats = F)

# Set plotting parameters
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

r <- rasters[[i]]

# Plot correlation raster with color scale
plot(r %>% terra::mask(., Hds_county_outline), 
     range = rn, 
     col = contrast.color(16),
     plg = list(title = "Spearman's Rank Correlation"),
     mar = c(0.5, 0.5, 0.1, 1), axes = FALSE, border = TRUE, alpha = 1, legend = FALSE)

# Overlay hillshade for terrain context
plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)

# Overlay county boundaries
plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)

# Overlay river layer
plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 1.5)

par(xpd = NA)  # allow drawing outside the plot region
# Add panel letter 'a' in the top-left corner (custom function `panel.let()`)
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = pic_order[i], font = 1, cex = 2)

# Add transparent-background legend for rivers and faults
legend(280000, 2980000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 1, bg = "transparent")

# Plot color legend manually (replacement for custom cscl() function)

# Define break points for color edges (ensure exact alignment with color bins)
breaks <- seq(rn[1], rn[2], length.out = 17)  # 16 color bins = 17 breaks

fields::image.plot(legend.only = TRUE,
                   zlim = rn,
                   col = contrast.color(16),
                   breaks = breaks,  # enforce correct bin edges
                   smallplot = c(0.87, 0.9, 0.2, 0.9),  # legend position
                   legend.lab = "Spearman's Rank Correlation",
                   legend.cex = 0.8,
                   axis.args = list(
                     at = lab,                     # tick positions
                     labels = names(lab),          # display as "-1", "-0.5", ..., "1"
                     cex.axis = 0.8
                   ),
                   legend.line = -1.8)
par(xpd = FALSE)

# Close the graphic device
dev.off()

# Clean up memory
gc()

}

# Dimensions
width_cm <- 7.95
height_cm <- 6.9
dpi <- 300
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

# Read in each plot
fig1 <- image_read_pdf("output_SI/sf10_beta_window_12km.pdf", density = dpi)
fig2 <- image_read_pdf("output_SI/sf10_beta_window_18km.pdf", density = dpi)
fig3 <- image_read_pdf("output_SI/sf10_beta_window_30km.pdf", density = dpi)
fig4 <- image_read_pdf("output_SI/sf10_beta_window_51km.pdf", density = dpi)

# Resize to match
fig1 <- image_resize(fig1, paste0(target_width, "x", target_height, "!"))
fig2 <- image_resize(fig2, paste0(target_width, "x", target_height, "!"))
fig3 <- image_resize(fig3, paste0(target_width, "x", target_height, "!"))
fig4 <- image_resize(fig4, paste0(target_width, "x", target_height, "!"))

# Combine into 2×2 layout
row1 <- image_append(c(fig1, fig2))
row2 <- image_append(c(fig3, fig4))
combined <- image_append(c(row1, row2), stack = TRUE)

# Save combined image
image_write(combined, path = "output_SI/sf10_beta_window_combined.png", format = "png")
image_write(combined, path = "output_SI/sf10_beta_window_combined.pdf", format = "pdf")

### =========================================================================
### Plot Supplementary Figure 11 - Residuals of a linear model of endemic 
### richness along elevation belts
### =========================================================================

# Load libraries for advanced plotting
library(plot3D)
library(lattice)         # for xyplot()
library(latticeExtra)    # for doubleYScale()
library(vip)
library(grid)
library(gridExtra)       # for grid.arrange() and viewports

# Function to create a double-y axis plot showing correlation between 
# a predictor (e.g., cost distance) and the residuals from SAR model
plot_5km <- function(i, summ) {
  
  # Fit a linear model between residuals and the selected variable `i`
  lm_5km <- lm(eval(parse(text = paste("resi~", i, sep = ""))), 
               data = summ %>%
                 scale() %>%           # scale all columns (standardize)
                 as.data.frame())      # convert back to data.frame
  
  # Extract model coefficient table and R-squared value
  lm_coeff <- summary(lm_5km)$coefficient
  lm_r <- summary(lm_5km)$r.squared
  lm_p <- round(lm_coeff[2, 4], 3)
  lm_est <- round(lm_coeff[2, 1], 3)
  
  # Format all text as plotmath using bquote with nested atop
  p_text <- if (lm_p == 0) {
    bquote(
      atop(
        atop("Coef = " * .(lm_est),
             R^2 == .(round(lm_r, 3))),
        italic(P) < 0.01
      )
    )
  } else {
    bquote(
      atop(
        atop("Coef = " * .(lm_est),
             R^2 == .(round(lm_r, 3))),
        italic(P) == .(lm_p)
      )
    )
  }
  
  # First Y-axis plot: residuals ~ elevation
  obj1 <- xyplot(resi ~ ele, data = summ, type = "l", lwd = 2, col = "steelblue",
                 scales = list(tck = 0), xlab = "Elevation (m)",
                 ylab = "Residuals",
                 panel = function(x, y, ...) {
                   panel.xyplot(x, y, ...);
                   # Add model annotation inside plot
                   ltext((min(summ$ele) + 800), 0.65 * max(summ$div), 
                         labels = p_text, cex = 1)
                 })
  
  # Second Y-axis plot: selected cost-distance or PC variable ~ elevation
  obj2 <- xyplot(eval(parse(text = paste(i, " ~ ele", sep = ""))), 
                 data = summ, type = "l", lwd = 2, col = "orange",
                 scales = list(tck = 0), xlab = "Elevation (m)",
                 # Dynamic y-label: handles Q vars or PC
                 ylab = ifelse(i == "PC",
                               i,
                               parse(text = paste("Cost~distance[", 
                                                  gsub("q_", "Q[", i), 
                                                  "]]", sep = "") %>%
                                       gsub(pattern = "n]]", replacement = "n]"))))
  
  # Combine the two plots using latticeExtra::doubleYScale()
  fig <- doubleYScale(obj1, obj2, add.ylab2 = TRUE, use.style = FALSE)
  
  return(fig)
}

# Load elevational summary data (rows = elev bands, cols = variables)
summ <- read.csv("./File/ele_summ.csv", row.names = 1)

# Fit a species-area relationship (SAR) model: log-log form
sar <- lm(log(div) ~ log(area), data = summ)

# Calculate residuals from SAR, but back-transform fitted values (exp)
summ$resi <- summ$div - exp(sar$fitted.values)

# Fit a second model: endemic richness vs all richness
sar <- lm(div ~ all_div, data = summ)
summ$resi_all <- residuals(sar)

# Rename columns for clarity
colnames(summ)[colnames(summ) == "q_50"] <- "Median"
colnames(summ)[colnames(summ) == "mean"] <- "Mean"

# Generate plots for each quantile of cost-distance
fig_1 <- plot_5km("Mean", summ)
fig_2 <- plot_5km("q_05", summ)
fig_3 <- plot_5km("q_25", summ)
fig_4 <- plot_5km("Median", summ)   # this one is plotted in the final output
fig_5 <- plot_5km("q_75", summ)
fig_6 <- plot_5km("q_95", summ)


pdf("output_SI/sf11.pdf", width = 17 / 2.54, height = 14 / 2.54, useDingbats = FALSE, pointsize = 7.5)
grid.arrange(fig_1,fig_2,fig_3,fig_4,fig_5,fig_6, ncol = 3,
             vp=viewport(width = 0.9, height=1, just = 0.55, gp = gpar(fontsize = 6, lwd = 0.8, cex = 0.5)))
KeyA<-list(space="right",
           lines=list(lwd = 2, col = c("steelblue", "orange")),
           text=list(c("Residuals","Cost distance")), padding.text=10)
draw.key(KeyA, draw = TRUE, vp = viewport(.9, .50), gp = gpar(fontsize = 3, lwd = 0.5,  cex = 0.2))
dev.off()

### =========================================================================
### Plot Supplementary Figure 12 - Endemic β-Diversity Along an Elevational Gradient 
### =========================================================================


### part B - Relative importance of predictors for beta diversity

# Load required raster layers for NMDS axes and beta diversity
NMDS1_map <- raster("./File_SI/NMDS1_map.tif")
NMDS2_map <- raster("./File_SI/NMDS2_map.tif")
NMDS3_map <- raster("./File/NMDS3_map.tif")
beta_map <- raster("./File_SI/beta_map.tif")
names(beta_map) <- "beta"

# Load DEM and project it to match beta raster
ele <- raster("./File/1000m_UTM47_DEM.tif")
ele <- ele %>% projectRaster(to = beta_map)

# Stack and convert rasters to a data.frame for analysis
beta_all <- stack(ele, NMDS1_map, NMDS2_map, NMDS3_map, beta_map) %>%
  rasterToPoints() %>%
  as.data.frame() %>%
  na.omit()
colnames(beta_all) <- c("x", "y", "ele", "NMDS1", "NMDS2", "NMDS3", "beta")

# Load elevational summaries
summ <- read.csv("./File/ele_summ.csv", row.names = 1)

# Prepare empty data frame for beta summaries by elevational band
beta_summ <- data.frame(
  ele = NA,
  beta = NA,
  NMDS1 = NA,
  NMDS2 = NA,
  NMDS3 = NA,
  Max = NA,
  Min = NA
)

# Loop over elevational bands and calculate mean beta and NMDS SDs
for (i in 10:49 * 100) {
  beta_ele <- beta_all %>%
    filter(ele > i & ele < (i + 100))
  
  beta_summ <- rbind(beta_summ,
                     c(i * 1,
                       mean(beta_ele$beta),
                       sd(beta_ele$NMDS1),
                       sd(beta_ele$NMDS2),
                       sd(beta_ele$NMDS3),
                       (mean(beta_ele$beta) + (sd(beta_ele$beta) / sqrt(length(beta_ele$beta)))),  # Max (mean + SE)
                       (mean(beta_ele$beta) - (sd(beta_ele$beta) / sqrt(length(beta_ele$beta))))   # Min (mean - SE)
                     ))
}

# Remove NA rows
beta_summ <- na.omit(beta_summ)

# Merge beta summary values into elevational summary
summ <- merge(summ, beta_summ, by = "ele")

# Variables to test as predictors
exp_var <- c("mean", "q_05", "q_25", "q_50", "q_75", "q_95")

# Compare R² for each cost-distance predictor variable
conn_R2 <- data.frame(
  exp = exp_var,
  R2 = NA
)
for (i in 1:length(exp_var)) {
  conn_lm <- lm(beta ~ ., data = summ[, c("beta", exp_var[i])]) %>% summary()
  conn_R2[i, 2] <- conn_lm$r.squared
}
# Identify the best cost-distance predictor
conn_ind <- conn_R2[which.max(conn_R2$R2), 1]

# Fit a multiple regression using best cost-distance predictor and environmental variables
lm_summ <- lm(beta ~ tem + I(tem^2) + velo + I(velo^2) + q_50, 
              data = summ[, c("beta", "q_50", "tem", "swb", "velo")] %>%
                scale() %>%
                as.data.frame())

# Extract adjusted R²
R2 <- summary(lm_summ)$adj.r.squared

# Compute relative importance using relaimpo package
imp_summ <- lm_summ %>% relaimpo::calc.relimp(rela = T)
imp <- imp_summ@lmg %>% sort(decreasing = TRUE)

# Format into grouped (linear/quadratic) dataframe
imp <- data.frame(
  Linear = imp[c("tem", "velo", "q_50")],
  Quadratic = imp[c("I(tem^2)", "I(velo^2)", "")]
)

# Reorder rows for plotting
imp <- imp[c(1, 3, 2), ]

# Extract coefficients corresponding to selected predictors (not used in plot here)
coef <- summary(lm_summ)$coefficients[, 1]
coef <- coef[names(imp)]

# Start PDF device for plotting
pdf("output_SI/sf12b.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

# Basic plotting setup
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Create barplot of relative importance
barplot(imp %>% as.matrix %>% t(), 
        col = c("steelblue1", "lightpink"), 
        names = imp %>% rownames() %>%
          gsub(pattern = "tem", replacement = "MTWQ") %>%
          gsub(pattern = "q_50", replacement = "Cost ~ Distance[Median]") %>%
          gsub(pattern = "velo", replacement = "ClimVelo") %>%
          parse(file = NULL, n = NULL),
        xlab = "Predictors",
        ylab = "Relavtive importance",
        legend = imp %>% colnames(),
        args.legend = list(bty = "n", x = "topright", inset = c(0, -0.2)))

# Add R² text annotation
text(
  x = nrow(imp) + 0.1, 0.5,
  labels = bquote(R^2 == .(round(summary(lm_summ)$r.squared, 3)))
)

# Add panel label 'b'
par(xpd = NA)
text(x = par("usr")[1] + -0.10 * diff(par("usr")[1:2]),
     y = par("usr")[3] + 1.15 * diff(par("usr")[3:4]),
     labels = "b", font = 1, cex = 2)
par(xpd = FALSE)

# Close device
dev.off()



### Part B - Elevational trend in beta diversity and cost distance

# Load tactile package for confidence interval plotting
library(tactile)

# Rename 'q_50' column to more descriptive name
colnames(summ)[colnames(summ) == "q_50"] <- "Median"

# Set target variable to be used for plotting
i <- "Median"

# Optional: Custom panel for confidence bands (not used here but defined)
my.panel.bands <- function(x, y, upper, lower, fill, col,
                           subscripts, ..., font, fontface)
{
  upper <- upper[subscripts]
  lower <- lower[subscripts]
  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                col = fill, border = FALSE,
                ...)
}

# Fit linear model between beta diversity and cost-distance (Median)
lm_5km <- lm(eval(parse(text = paste("beta~", i))), 
             data = summ %>%
               scale() %>%           # standardize variables
               as.data.frame())      # convert back to data.frame

# Extract model summary statistics
lm_coeff <- summary(lm_5km)$coefficient
lm_r <- summary(lm_5km)$r.squared

# Format coefficient and p-value
lm_est <- round(lm_coeff[2, 1], 3)
lm_r <- round(lm_r, 3)
lm_pval <- round(lm_coeff[2, 4], 3)

# Format display text for model stats (as expression)
p_text <- if (lm_pval == 0) {
  bquote(
    atop(
      atop("Coef =" ~ .(lm_est),
           R^2 == .(lm_r)),
      italic(P) < 0.01
    )
  )
} else {
  bquote(
    atop(
      atop("Coef =" ~ .(lm_est),
           R^2 == .(lm_r)),
      italic(P) == .(lm_pval)
    )
  )
}

# Load latticeExtra for double axis plotting and ltext
library(latticeExtra)

# Panel 1: Plot beta diversity with shaded confidence interval
obj1 <- xyplot(beta ~ ele, data = summ, type = "l", lwd = 2, col = "steelblue",
               upper = summ$Max, lower = summ$Min,
               scales = list(tck = 0), xlab = "Elevation (m)",
               ylab = "Beta diversity",
               panel = function(x, y, ...) {
                 panel.ci(x = x, y = y, upper = summ$Max, lower = summ$Min, 
                          alpha = 0.15, grid = FALSE)
                 panel.xyplot(x, y, ...)
                 ltext((min(summ$ele) + 800), 0.95 * max(summ$beta), 
                       labels = p_text, cex = 1)
               })

# Panel 2: Plot cost-distance (Median) on secondary y-axis
obj2 <- xyplot(eval(parse(text = paste(i, " ~ ele",  sep = ""))), 
               data = summ, type = "l", lwd=2, col="orange",
               scales = list(tck = 0), xlab = "Elevation (m)",
               ylab = ifelse(i == "PC",
                             i,
                             parse(text = paste("Cost~Distance[", gsub("q_", "Q[", i) , "]]" , sep = "")%>%
                                     gsub(pattern = "n]]", replacement = "n]"))))

# Combine both panels using dual y-axis
fig <- doubleYScale(obj1, obj2, add.ylab2 = TRUE, use.style = FALSE)

# Start PDF device for publication figure
pdf("output_SI/sf12a.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

# Display the composite plot with layout adjustments
grid.arrange(fig, ncol = 1,
             vp = viewport(width = 1, height = 1, just = 0.5, 
                           gp = gpar(fontsize = 7, lwd = 1, cex = 0.5)))

# Legend for the two lines (beta diversity and cost distance)
KeyA <- list(space = "right",
             lines = list(lwd = 2, col = c("steelblue", "orange")),
             text = list(c("Beta diversity\n(w/SE)", "Cost distance")), 
             padding.text = 5)

# Panel label 'a'
KeyB <- list(space = "right",
             text = list('a'), 
             padding.text = 5)

# Draw both legends using grid viewports
draw.key(KeyA, draw = TRUE, vp = viewport(.75, .7, gp = gpar(fontsize = 5, lwd = 1)))
draw.key(KeyB, draw = TRUE, vp = viewport(0.05, .95, gp = gpar(fontsize = 15, lwd = 1)))

# Close plotting device
dev.off()

### Figure combine

# Define file paths to the two PDF figures
fig1_path <- "output_SI/sf12a.pdf"
fig2_path <- "output_SI/sf12b.pdf"

# Define desired physical size (in cm) and resolution (in DPI)
width_cm <- 7.95
height_cm <- 6.9
dpi <- 300

# Convert dimensions to pixels
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

# Read the first page of each PDF as an image
fig1 <- image_read_pdf(fig1_path, density = dpi) %>%
  image_resize(paste0(target_width, "x", target_height, "!"))

fig2 <- image_read_pdf(fig2_path, density = dpi) %>%
  image_resize(paste0(target_width, "x", target_height, "!"))

# Combine the two images side-by-side (horizontal layout)
combined_image <- image_append(c(fig1, fig2))

# Save the combined image to a file
image_write(combined_image, path = "output_SI/sf12_combined.png", format = "png")
image_write(combined_image, path = "output_SI/sf12_combined.pdf", format = "pdf")

### =========================================================================
### Plot Supplementary Figure 13 - NDMS axis 1 and 2 
### =========================================================================

# change raster format
NMDS1_map <- rast(NMDS1_map)
NMDS2_map <- rast(NMDS2_map)

# Scale the raster values to range [-1, 1] for color normalization
scaled_raster1 <- ((NMDS1_map - minmax(NMDS1_map)[1]) / 
                    (minmax(NMDS1_map)[2] - minmax(NMDS1_map)[1]) - 0.5) * 2

scaled_raster2 <- ((NMDS2_map - minmax(NMDS2_map)[1]) / 
                     (minmax(NMDS2_map)[2] - minmax(NMDS2_map)[1]) - 0.5) * 2

# Define display range for the scaled values
rn = c(-1, 1)

# Define tick positions and corresponding labels for legend
lab = c('-1' = -1,
        "-0.5" = -0.5,
        "0" = 0,
        "0.5" = 0.5,
        "1" = 1)

# Create a diverging color palette for the NDMS axis
NDMS_col <- colorRampPalette(c(
  '#543005',
  '#8c510a',
  '#bf812d',
  '#dfc27d',
  '#f6e8c3',
  '#f5f5f5',
  '#c7eae5',
  '#80cdc1',
  '#35978f',
  '#01665e',
  '#003c30'
))

### part A axis 1

# Save to publication-quality PDF
pdf("output_SI/sf13a.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, pointsize = 7.5, useDingbats = FALSE)

# Set plot layout and base font size
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Plot the scaled NDMS raster with no legend
plot(scaled_raster1 %>% terra::mask(., Hds_county_outline), 
     range = rn, 
     col = NDMS_col(17),  # 17 color levels
     plg = list(title = "NDMS1"),  # internal label, not used here
     mar = c(0.5, 0.5, 0.1, 1), axes = FALSE, border = TRUE, alpha = 0.8, legend = FALSE)

# Add hillshade background for terrain effect
plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)

# Overlay county boundary
plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)

# Overlay river vector in ETH blue
plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 2)

# Overlay fault lines as dashed black lines
plot(faults, col = 'black', add = TRUE, border = FALSE, lty = 2, lwd = 2)

# Allow drawing outside plot margin for annotations
par(xpd = NA)

# Add panel letter label "d" in upper left corner
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = "a", font = 1, cex = 2)

# Add map legend for river and fault overlays
legend(280000, 2980000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 1, bg = "transparent")
legend(280000, 2930000, legend = "Faults", col = "black", lty = 2, lwd = 2, box.lty = 0, cex = 1, bg = "transparent")

# Create breaks to align color edges with the data scale
breaks <- seq(rn[1], rn[2], length.out = 17)  # 16 color breaks = 17 edges

# Plot vertical color scale legend
fields::image.plot(legend.only = TRUE,
                   zlim = rn,
                   col = NDMS_col(16),        # pass one fewer color than breaks
                   breaks = breaks,           # ensures precise alignment of color blocks
                   smallplot = c(0.87, 0.9, 0.2, 0.9),  # controls legend size and placement
                   legend.lab = "NDMS Axis 1 Values",  # legend title
                   legend.cex = 0.8,
                   axis.args = list(
                     at = lab,                # tick mark positions
                     labels = names(lab),     # corresponding labels
                     cex.axis = 0.8
                   ),
                   legend.line = -1.8)  # position of legend label relative to axis

# Turn off plotting outside region
par(xpd = FALSE)

# Close graphics device
dev.off()

# Garbage collection to clear memory
gc()


### part B axis 2 

pdf("output_SI/sf13b.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, pointsize = 7.5, useDingbats = FALSE)

# Set plot layout and base font size
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Plot the scaled NDMS raster with no legend
plot(scaled_raster2 %>% terra::mask(., Hds_county_outline), 
     range = rn, 
     col = NDMS_col(17),  # 17 color levels
     plg = list(title = "NDMS2"),  # internal label, not used here
     mar = c(0.5, 0.5, 0.1, 1), axes = FALSE, border = TRUE, alpha = 0.8, legend = FALSE)

# Add hillshade background for terrain effect
plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)

# Overlay county boundary
plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)

# Overlay river vector in ETH blue
plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 2)

# Overlay fault lines as dashed black lines
plot(faults, col = 'black', add = TRUE, border = FALSE, lty = 2, lwd = 2)

# Allow drawing outside plot margin for annotations
par(xpd = NA)

# Add panel letter label "d" in upper left corner
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = "b", font = 1, cex = 2)

# Add map legend for river and fault overlays
legend(280000, 2980000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 1, bg = "transparent")
legend(280000, 2930000, legend = "Faults", col = "black", lty = 2, lwd = 2, box.lty = 0, cex = 1, bg = "transparent")

# Create breaks to align color edges with the data scale
breaks <- seq(rn[1], rn[2], length.out = 17)  # 16 color breaks = 17 edges

# Plot vertical color scale legend
fields::image.plot(legend.only = TRUE,
                   zlim = rn,
                   col = NDMS_col(16),        # pass one fewer color than breaks
                   breaks = breaks,           # ensures precise alignment of color blocks
                   smallplot = c(0.87, 0.9, 0.2, 0.9),  # controls legend size and placement
                   legend.lab = "NDMS Axis 2 Values",  # legend title
                   legend.cex = 0.8,
                   axis.args = list(
                     at = lab,                # tick mark positions
                     labels = names(lab),     # corresponding labels
                     cex.axis = 0.8
                   ),
                   legend.line = -1.8)  # position of legend label relative to axis

# Turn off plotting outside region
par(xpd = FALSE)

# Close graphics device
dev.off()

# Garbage collection to clear memory
gc()

### combine

# Define file paths to the two PDF figures
fig1_path <- "output_SI/sf13a.pdf"
fig2_path <- "output_SI/sf13b.pdf"

# Define desired physical size (in cm) and resolution (in DPI)
width_cm <- 7.95
height_cm <- 6.9
dpi <- 300

# Convert dimensions to pixels
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

# Read the first page of each PDF as an image
fig1 <- image_read_pdf(fig1_path, density = dpi) %>%
  image_resize(paste0(target_width, "x", target_height, "!"))

fig2 <- image_read_pdf(fig2_path, density = dpi) %>%
  image_resize(paste0(target_width, "x", target_height, "!"))

# Combine the two images side-by-side (horizontal layout)
combined_image <- image_append(c(fig1, fig2))

# Save the combined image to a file
image_write(combined_image, path = "output_SI/sf13_combined.png", format = "png")
image_write(combined_image, path = "output_SI/sf13_combined.pdf", format = "pdf")

### =========================================================================
### Supplementary Figure 14 - elevational predictor importance for endemic
### richness and beta diversity
### =========================================================================

# === Load input rasters ===
ele <- raster("./File/1000m_UTM47_DEM.tif")
end_div <- raster("./File/end_div.tif")
Flux <- raster("./File/flux_mean_d10.tif")
velo <- raster("./File/Velo.tif")
tem <- raster("./File/tem.tif")
swb <- raster("./File/swb.tif")

### === Part A: Endemic diversity ===
# === Extract values and construct dataframe for modeling ===
flux_all <- data.frame(
  ele = values(ele),
  end_div = values(end_div),
  tem = values(tem),
  swb = values(swb),
  Flux = values(Flux),
  velo = abs(values(velo))  # use absolute values for climate velocity
) %>%
  cbind(coordinates(ele)) %>%  # add spatial coordinates
  na.omit() %>%
  filter(is.finite(Flux))  # remove rows with non-finite Flux values

# === Fit Poisson GLM model for endemic richness ===
all_glm <- glm(
  end_div ~ tem + swb + log(velo) + Flux + I(tem^2) + I(swb^2) + I(log(velo)^2) + I(Flux^2),
  family = "poisson",
  data = flux_all
)

# === Calculate variable importance and format labels ===
all_imp <- caret::varImp(all_glm)
all_imp <- structure(all_imp[, 1], names = rownames(all_imp))
all_imp <- all_imp / sum(all_imp)
all_imp <- sort(all_imp, decreasing = TRUE)

# Rename terms for clarity in legend
names(all_imp)[names(all_imp) == "log(velo)"] <- "ClimVelo"
names(all_imp)[names(all_imp) == "I(log(velo)^2)"] <- "I(ClimVelo^2)"

# Separate linear and quadratic terms
all_imp <- data.frame(
  Linear = all_imp[c("tem", "Flux", "ClimVelo", "swb")],
  Quadratic = all_imp[c("I(tem^2)", "I(Flux^2)", "I(ClimVelo^2)", "I(swb^2)")]
)

# === Plot output for richness ===
pdf("output_SI/sf14a.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

par(mfrow = c(1,1), oma = c(0,0,0,0), ps = 8, cex = 1)

barplot(t(as.matrix(all_imp)), 
        col = c("steelblue1", "lightpink"),
        xlab = "Predictors",
        ylab = "Relative importance",
        names = rownames(all_imp) %>%
          gsub("tem", "MTWQ", .) %>%
          gsub("swb", "SWB", .) %>%
          parse(text = .),
        legend = colnames(all_imp),
        args.legend = list(bty = "n", x = 5, y = 0.6))

# Add model explanatory power as D²
text(
  x = nrow(all_imp) + 0.1,
  y = 0.68 * max(rowSums(all_imp)),
  labels = bquote(D^2 == .(round(ecospat::ecospat.adj.D2.glm(all_glm), 3))),
  cex = 1
)

title(main = "Endemic species α-diversity")

# Panel label
par(xpd = NA)
text(x = par("usr")[1] + -0.10 * diff(par("usr")[1:2]),
     y = par("usr")[3] + 1.15 * diff(par("usr")[3:4]),
     labels = "a", font = 1, cex = 2)
par(xpd = FALSE)

dev.off()

### === Part B: Beta diversity ===

# === Load and align predictors to beta diversity map extent ===
beta_map <- raster("./File_SI/beta_map.tif")
ele1 <- resample(ele, beta_map) %>% rast() %>% raster::mask(Hds_county_outline) %>% raster()
tem1 <- resample(tem, beta_map) %>% rast() %>% raster::mask(Hds_county_outline) %>% raster()
swb1 <- resample(swb, beta_map) %>% rast() %>% raster::mask(Hds_county_outline) %>% raster()
Flux1 <- resample(Flux, beta_map) %>% rast() %>% raster::mask(Hds_county_outline) %>% raster()
velo1 <- resample(velo, beta_map) %>% rast() %>% raster::mask(Hds_county_outline) %>% raster()

# === Combine into dataframe ===
flux_all <- data.frame(
  ele = values(ele1),
  beta = values(beta_map),
  tem = values(tem1),
  swb = values(swb1),
  Flux = values(Flux1),
  velo = abs(values(velo1))
) %>%
  cbind(coordinates(ele1)) %>%
  na.omit() %>%
  filter(is.finite(Flux))

# === Fit binomial GLM for beta diversity ===
all_glm <- glm(
  beta ~ tem + swb + log(velo) + Flux + I(tem^2) + I(swb^2) + I(log(velo)^2) + I(Flux^2),
  family = "binomial",
  data = flux_all
)

# === Extract variable importance ===
all_imp <- caret::varImp(all_glm)
all_imp <- structure(all_imp[,1], names = rownames(all_imp))
all_imp <- all_imp / sum(all_imp)
all_imp <- sort(all_imp, decreasing = TRUE)
names(all_imp)[names(all_imp) == "log(velo)"] <- "ClimVelo"
names(all_imp)[names(all_imp) == "I(log(velo)^2)"] <- "I(ClimVelo^2)"

# Separate and order importance values
all_imp <- data.frame(
  Linear = all_imp[c("tem", "Flux", "ClimVelo", "swb")],
  Quadratic = all_imp[c("I(tem^2)", "I(Flux^2)", "I(ClimVelo^2)", "I(swb^2)")]
)
all_imp <- all_imp[c(1, 3, 2, 4), ]  # reordering for consistent bar order

# === Plot output for beta diversity ===
pdf("output_SI/sf14b.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

par(mfrow = c(1,1), oma = c(0,0,0,0), ps = 8, cex = 1)

barplot(t(as.matrix(all_imp)), 
        col = c("steelblue1", "lightpink"),
        xlab = "Predictors",
        ylab = "Relative importance",
        names = rownames(all_imp) %>%
          gsub("tem", "MTWQ", .) %>%
          gsub("swb", "SWB", .) %>%
          parse(text = .),
        legend = colnames(all_imp),
        args.legend = list(bty = "n", x = 5, y = 0.7))

# Add adjusted D² for model fit
text(
  x = nrow(all_imp) + 0.1,
  y = 0.68 * max(rowSums(all_imp)),
  labels = bquote(D^2 == .(round(ecospat::ecospat.adj.D2.glm(all_glm), 3))),
  cex = 1
)

title(main = "Endemic species β-diversity")

# Add subfigure label
par(xpd = NA)
text(x = par("usr")[1] + -0.10 * diff(par("usr")[1:2]),
     y = par("usr")[3] + 1.15 * diff(par("usr")[3:4]),
     labels = "b", font = 1, cex = 2)
par(xpd = FALSE)

dev.off()

### combine

# Define file paths to the two PDF figures
fig1_path <- "output_SI/sf14a.pdf"
fig2_path <- "output_SI/sf14b.pdf"

# Define desired physical size (in cm) and resolution (in DPI)
width_cm <- 7.95
height_cm <- 6.9
dpi <- 300

# Convert dimensions to pixels
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

# Read the first page of each PDF as an image
fig1 <- image_read_pdf(fig1_path, density = dpi) %>%
  image_resize(paste0(target_width, "x", target_height, "!"))

fig2 <- image_read_pdf(fig2_path, density = dpi) %>%
  image_resize(paste0(target_width, "x", target_height, "!"))

# Combine the two images side-by-side (horizontal layout)
combined_image <- image_append(c(fig1, fig2))

# Save the combined image to a file
image_write(combined_image, path = "output_SI/sf14_combined.png", format = "png")
image_write(combined_image, path = "output_SI/sf14_combined.pdf", format = "pdf")