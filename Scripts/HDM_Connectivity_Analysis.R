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
### Plot endemic richness Figure 1 A
### =========================================================================

# Define legend tick positions and labels
rn <- c(0, 773)
lab_vals <- c(0, 150, 300, 450, 600, 750)
lab_names <- as.character(lab_vals)
lab_vals <- lab_vals[lab_vals >= rn[1] & lab_vals <= rn[2]]

# Start PNG device
#png("output/F1A_endemic_richness.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)
pdf("output/F1A_endemic_richness.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)  # for publication-ready PDF

par(mar = c(0.1, 0.1, 0.1, 4.2), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Plot richness raster without axis/legend
plot(end_div,col =col_vec,legend=F,axes=F,mar=c(0.1,0.1,0.1,4.2), border=F)

# Overlay hillshade to add terrain perception
plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)

# Overlay river vector
plot(river, col = ETH_blue_old, lwd = 1.5, add = TRUE)

# Overlay county outlines
plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)

# Define color breaks matching col_vec
n_col <- length(col_vec)
breaks <- seq(rn[1], rn[2], length.out = n_col + 1)

# Add vertical color legend with aligned ticks
image.plot(legend.only = TRUE,
           zlim = rn,
           col = col_vec,
           breaks = breaks,
           smallplot = c(0.87, 0.9, 0.2, 0.9),  # position legend at right edge
           legend.lab = "Endemic richness",
           axis.args = list(at = lab_vals, labels = lab_names, cex.axis = 0.8),
           legend.line = -2)

# Add subpanel label 'a'
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = "a", font = 1, cex = 2)

# Add legend for river
legend(-150000, 2900000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 0.9)

dev.off()
gc()

### =========================================================================
### Plot elevation profile Figure 1 B
### =========================================================================

# Load elevational summary data
summ <- read.csv("./File/ele_summ.csv", row.names = 1)

# Extract selected columns for plotting (area, endemic, overall richness)
summ.1c <- summ[, c(8,9,15)]

# Define min-max scaling function (normalize values between 0 and 1)
min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Apply scaling
summ.1c <- as.data.frame(lapply(summ.1c, min_max_scale))
summ.1c$ele <- summ$ele  # re-attach elevation column for x-axis

# Optional: proportion endemic
# summ.1c$por <- summ$div / summ$all_div

# Prepare bar data for area
summ.1c.bar <- summ.1c$area %>% as.matrix()
new_row_names <- summ.1c$ele
rownames(summ.1c.bar) <- new_row_names

# Output image
#png("output/f1b.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)

pdf("output/f1b.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

par(mfrow = c(1,1), oma = c(0,0,0,0), ps = 8, cex = 1)

# Create bar and line plot using ggplot2
ggplot(summ.1c, aes(x = ele)) +
  geom_bar(aes(y = area, fill = 'Area'), stat = "identity", color = "white") +
  geom_line(aes(y = div, color = "Endemic"), size = 1) +
  geom_line(aes(y = all_div, color = "Overall"), size = 1) +
  #geom_line(aes(y = por, color = "%endemic"), size = 1) +
  xlab("Elevation (m)") +
  ylab("Standardized value") +
  labs(color = "Richness type", fill = 'Area' , tag = "b") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.3, 'cm'), 
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'), 
        legend.box.margin = margin(-10, -10, -10, -10)) +
  scale_color_manual(values = c("Endemic" = "#0099CC", "Overall" = "#FF9933", 'Area' = 'lightpink'))

dev.off()

### =========================================================================
### Plot glaciation Figure 1 C
### =========================================================================

# Define elevation legend range
rn <- c(0, 7500)
lab_vals <- c(0, 1500, 3000, 4500, 6000, 7500)
lab_names <- as.character(lab_vals)

# png("output/glaciation_f1c.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)
pdf("output/glaciation_f1c.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

par(mar = c(0.1, 0.1, 0.1, 4.2), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Plot raw elevation
plot(ele_raw,col =terrain.colors(50),legend=F,axes=F,mar=c(0.1,0.1,0.1,4.2), plg = list(title = "Elevation(m)"), border=T , alpha = 0.9)

# Add hillshade
plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)

# Add LGM glaciation extent
plot(glacial, col = ETH_blue_new, lwd = 1.5, add = TRUE, border = FALSE, legend=F)

# Add rivers
plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 1.5)

# Overlay county outlines
plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)

# Define color vector and breaks for elevation
elev_col <- terrain.colors(50)
n_col <- length(elev_col)
breaks <- seq(rn[1], rn[2], length.out = n_col + 1)

# Add elevation legend with properly aligned ticks
image.plot(legend.only = TRUE,
           zlim = rn,
           col = elev_col,
           breaks = breaks,
           smallplot = c(0.87, 0.9, 0.2, 0.9),
           legend.lab = "Elevation (m)",
           axis.args = list(at = lab_vals, labels = lab_names, cex.axis = 0.7),
           legend.line = -2)

# Panel label 'c'
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = "c", font = 1, cex = 2)

# Add legends for glacier and rivers
legend(-150000, 2900000, legend = 'LGM glacier', fill = ETH_blue_new, border = 'white', box.lty = 0, cex = 0.9)
legend(-150000, 2800000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 0.9)

dev.off()
gc()

### =========================================================================
### Plot elevation classes Figure 1 D
### =========================================================================

# Reclassify elevation into 5 bands
ele_raw1 <- ele_raw
ele_raw1[ele_raw1 < 1800] <- 1
ele_raw1[ele_raw1 >= 1800 & ele_raw1 < 2800] <- 2
ele_raw1[ele_raw1 >= 2800 & ele_raw1 < 3800] <- 3
ele_raw1[ele_raw1 >= 3800 & ele_raw1 < 4600] <- 4
ele_raw1[ele_raw1 >= 4600] <- 5

# Define 5-level color palette
col_terrain2 <- colorRampPalette(c("#008837", "#91cf60", "yellow", "#fc8d59", '#392a48'))(5)

# png("output/ele_f1d.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)
pdf("output/ele_f1d.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

par(mar = c(0.1, 0.1, 0.1, 4.2), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Plot classified elevation
plot(ele_raw1,col =col_terrain2,legend=F,axes=F,mar=c(0.1,0.1,0.1,4.2), plg = list(title = "Elevation(m)"), border=T , alpha = 0.9)


# Add hillshade
plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)

# Add rivers
plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 1.5)

# Overlay county outlines
plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)

# Subfigure label 'd'
par(xpd = NA)
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = "d", font = 1, cex = 2)

# Add legend for rivers and elevation classes
legend(-150000, 2990000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 0.8, bg = NA)
legend(-150000, 2910000, legend = c('below 1800m', '1800–2800m', '2800–3800m', '3800–4600m', 'above 4600m'), 
       fill = col_terrain2, border = 'white', box.lty = 0, cex = 0.8, bg = NA)
par(xpd = FALSE)

dev.off()
gc()

### =========================================================================
### Figure 1 combination
### =========================================================================

# figure 1 combine
# Read in PDF pages and convert to images
fig_a <- image_read_pdf("output/F1A_endemic_richness.pdf", density = 300)
fig_b <- image_read_pdf("output/f1b.pdf", density = 300)
fig_c <- image_read_pdf("output/glaciation_f1c.pdf", density = 300)
fig_d <- image_read_pdf("output/ele_f1d.pdf", density = 300)

# Resize images to same dimensions (optional, good for alignment)
# Desired output size in cm
width_cm <- 7.95
height_cm <- 6.9

# Desired resolution (DPI)
dpi <- 300

# Convert to pixels
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

fig_a <- image_resize(fig_a, paste0(target_width, "x", target_height, "!"))
fig_b <- image_resize(fig_b, paste0(target_width, "x", target_height, "!"))
fig_c <- image_resize(fig_c, paste0(target_width, "x", target_height, "!"))
fig_d <- image_resize(fig_d, paste0(target_width, "x", target_height, "!"))

# Combine in 2 rows, 2 columns
row1 <- image_append(c(fig_a, fig_b))  # horizontal append: a | b
row2 <- image_append(c(fig_c, fig_d))  # horizontal append: c | d
final <- image_append(c(row1, row2), stack = TRUE)  # vertical stack

# Save as final figure
image_write(final, path = "output/Figure1_combined.png", format = "png", density = 300)
# Optional: also save as PDF
image_write(final, path = "output/Figure1_combined.pdf", format = "pdf")

### =========================================================================
### Plot flux vs richness Spearman figure 2 A
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
res_rel_imp <- rast('./Results/res_rel_imp_average_50.tif')

# Start PNG output
#png("output/richvsflux_f2a.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)
pdf("output/richvsflux_f2a.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

par(mar = c(0.1, 0.1, 0.1, 4.2), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Mask correlation raster with study boundary
res_masked <- mask(res_rel_imp, Hds_county_outline)

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
     labels = "a", font = 1, cex = 2)

# Add river legend
legend(-150000, 2800000,
       legend = "Rivers",
       col = ETH_blue_old,
       lty = 1, lwd = 2,
       box.lty = 0, cex = 0.9)

dev.off()
gc()

### =========================================================================
### Plot  figure 2 B correlation on elevation 
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

# Export selected plot (Median vs residuals) to PNG
# png("output/f2b.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7)
pdf("output/f2b.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

# Display the Median cost-distance plot using grid graphics
grid.arrange(fig_4, ncol = 1,
             vp = viewport(width = 1, height = 1, just = 0.5, 
                           gp = gpar(fontsize = 7, lwd = 1, cex = 0.5)))

# Create custom keys (legends) for the lines
KeyA <- list(space = "right",
             lines = list(lwd = 2, col = c("steelblue", "orange")),
             text = list(c("Residuals", "Cost distance")),
             padding.text = 5)

# Subfigure label key ('b')
KeyB <- list(space = "right",
             text = list('b'), 
             padding.text = 5)

# Draw the two keys in specific positions using viewports
draw.key(KeyA, draw = TRUE, vp = viewport(.25, .7, gp = gpar(fontsize = 5, lwd = 1)))
draw.key(KeyB, draw = TRUE, vp = viewport(0.05, .95, gp = gpar(fontsize = 15, lwd = 1)))

dev.off()
gc()

### =========================================================================
### Plot  figure 2 c linear model on elevation bands
### =========================================================================

# Load elevational summary data (e.g. for elevational bands)
summ <- read.csv("./File/ele_summ.csv", row.names = 1)

# Fit a species-area relationship (SAR) model using log-log form
sar <- lm(log(div) ~ log(area), data = summ)

# Back-transform fitted values and calculate SAR residuals
# This gives the difference between observed and expected richness
summ$resi <- summ$div - exp(sar$fitted.values)

# Fit a second model predicting endemic richness using total richness
sar <- lm(div ~ all_div, data = summ)
summ$resi_all <- residuals(sar)  # save residuals for alternative use

# Fit OLS model predicting SAR residuals using environmental predictors
all_glm <- lm(
  resi ~ q_50 + tem + I(tem^2) + velo + I(velo^2),
  data = summ
)  

# Compute relative importance of predictors using LMG metric (from relaimpo package)
all_imp <- all_glm %>% relaimpo::calc.relimp(rela = TRUE)
all_imp <- all_imp@lmg  # extract LMG contributions
all_imp <- all_imp / sum(all_imp)  # normalize to sum to 1
all_imp <- sort(all_imp, decreasing = TRUE)  # rank by importance

# Rename predictors for clarity in final plot
names(all_imp)[names(all_imp) == "velo"] <- "ClimVelo"
names(all_imp)[names(all_imp) == "I(velo^2)"] <- "I(ClimVelo^2)"

# Format for barplot: separate linear and quadratic terms
all_imp <- data.frame(
  Linear = all_imp[c("q_50", "tem", "ClimVelo")],
  Quadratic  = c(0, all_imp["I(tem^2)"], all_imp["I(ClimVelo^2)"])
)

# Output figure
#png("output/f2c.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7)
pdf("output/f2c.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, useDingbats = FALSE, pointsize = 7.5)

# Plot stacked barplot of relative importance
barplot(
  all_imp %>% as.matrix() %>% t(),  # transpose for stacking
  col = c("steelblue1", "lightpink"),  # colors for linear and quadratic
  names = all_imp %>% rownames() %>%  # rename x-axis labels using plotmath
    gsub(pattern = "tem", replacement = "MTWQ") %>%
    gsub(pattern = "swb", replacement = "SWB") %>%
    gsub(pattern = "q_50", replacement = "Cost ~ Distance[Median]") %>%
    gsub(pattern = "q_75", replacement = "Cost ~ Distance[Q[75]]") %>% 
    gsub(pattern = "q_25", replacement = "Cost ~ Distance[Q[25]]") %>%
    gsub(pattern = "q_05", replacement = "Cost ~ Distance[Q[05]]") %>% 
    gsub(pattern = "q_95", replacement = "Cost ~ Distance[Q[95]]") %>%
    gsub(pattern = "mean", replacement = "Cost ~ Distance[Mean]") %>%
    parse(file = NULL, n = NULL),  # parse into expressions for plotmath labels
  xlab = "Predictors",
  ylab = "Relavtive importance",  # note: typo in original — should be "Relative importance"
  legend = all_imp %>% colnames(),  # add legend labels
  args.legend = list(bty = "n", x = 3.6, y = 0.6)  # legend placement
)

# Add R-squared annotation to the plot
text(
  x = nrow(all_imp) + 0.1,
  y = 0.75 * max(all_imp),
  labels = bquote(R^2 == .(round(summary(all_glm)$r.squared, 3)))
)

# Add figure panel label (e.g. "c") with slight offset and larger text
par(xpd=NA)
text(x = par("usr")[1] + -0.10 * diff(par("usr")[1:2]),
     y = par("usr")[3] + 1.15 * diff(par("usr")[3:4]),
     labels = "c", font = 1, cex = 2)
par(xpd=F)

# Close device and clean up
dev.off()

### =========================================================================
### Figure 2 combination
### =========================================================================

# figure 2 combine
# Read in PDF pages and convert to images
# Read in the images
f2a <- image_read_pdf("output/richvsflux_f2a.pdf", density = 300)
f2b <- image_read_pdf("output/f2b.pdf", density = 300)
f2c <- image_read_pdf("output/f2c.pdf", density = 300)

# Resize images to same dimensions (optional, good for alignment)
# Desired output size in cm
width_cm <- 7.95
height_cm <- 6.9

# Desired resolution (DPI)
dpi <- 300

# Convert to pixels
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

f2a <- image_resize(f2a, paste0(target_width, "x", target_height, "!"))
f2b <- image_resize(f2b, paste0(target_width, "x", target_height, "!"))
f2c <- image_resize(f2c, paste0(target_width, "x", target_height, "!"))

# Stack vertically
figure2 <- image_append(c(f2a, f2b, f2c), stack = TRUE)

# Save final combined figure
image_write(figure2, path = "output/Figure2_combined.png", format = "png", density = 300)
# Optional: save as PDF
image_write(figure2, path = "output/Figure2_combined.pdf", format = "pdf")

### ========================================================================
### Data prep for figure 3
### ========================================================================

# Load beta-diversity raster classified into categorical composition types
beta_cata <- rast("./File/beta_cata.tif")

# Define extent of the hotspot region (bounding box for cropping or plotting)
hotspot <- extent(360000, 880000, 2850000, 3500000)

# Define a custom 12-color palette for the 12 composition categories
my_color4 <- colorRampPalette(colors = c(
  '#73ABAB', '#BBE6E4',
  '#2FB8DA', '#074373',
  '#F57600', '#C44E00',
  '#A36912', '#663A00',
  '#E0E000', '#E0BB00',
  '#ff0505', '#b51936'))

### ========================================================================
### figure 3a endemic species assemblage at hotspot area
### ========================================================================

# Save output to PNG
# png("output/beta_cata_fig3a.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)

# Or save to PDF for vector graphics
pdf("output/beta_cata_fig3a.pdf", width = 7.95/2.54, height = 6.9/2.54, pointsize = 7.5, useDingbats = FALSE)

# Set graphical parameters: single panel, no outer margins, font size, base cex
par(mfrow = c(1,1), oma = c(0,0,0,0), ps = 8, cex = 1)

# Plot classified beta-diversity map with 12 color categories
plot(beta_cata, 
     col = my_color4(12),                        # apply color scale
     plg = list(title = "composition catagories"), # (note: typo in 'catagories')
     mar = c(0.5, 0.5, 0.1, 1),                  # inner margins
     axes = FALSE, border = TRUE, alpha = 1, legend = FALSE)

# Add hillshade to give terrain context
 plot(hillsh, bty = "n", legend = FALSE, axes = FALSE, col = rev(hillcol(50)), maxcell = 1e8, add = TRUE)

# Overlay county outlines
plot(Hds_county_outline, add = TRUE, lwd = 1, color = NA)

# Overlay rivers in ETH blue
plot(river, col = ETH_blue_old, add = TRUE, border = FALSE, lwd = 1.5)

# Overlay geological faults in dashed black lines
plot(faults, col = 'black', add = TRUE, border = FALSE, lty = 2, lwd = 2)

# Enable plotting outside the plot region for panel letter
par(xpd = NA)

# Add panel letter 'a' in the top-left corner (custom function `panel.let()`)
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
     labels = "a", font = 1, cex = 2)

# Disable outside plotting again
par(xpd = FALSE)

# Add transparent-background legend for rivers and faults
par(xpd = NA)
legend(280000, 2980000, legend = "Rivers", col = ETH_blue_old, lty = 1, lwd = 2, box.lty = 0, cex = 1, bg = "transparent")
legend(280000, 2930000, legend = "Faults", col = "black", lty = 2, lwd = 2, box.lty = 0, cex = 1, bg = "transparent")
par(xpd = FALSE)

# Close the graphical device
dev.off()

# Run garbage collection to clear memory
gc()

### ========================================================================
### categories fig3a-legend
### ========================================================================

######################
library(dendextend)
#####################
tree1 <- readRDS("./File/tree1.rds")

ForceUltrametric <- function(n) {
  if (is.leaf(n)) {
    # if object is a leaf, adjust height attribute
    attr(n, "height") <- 0.35
  }
  return(n)
}

tree1 <- dendrapply(X = tree1,
                    FUN = function(x) ForceUltrametric(x))

#write.tree(hc_phylo0, file = "12_cata_slice_tree.newick")

labels(tree1) <- as.character(seq(1:12))

pdf("output/beta_cata_fig3ab.pdf", width = 3.95/2.54, height = 6.9/2.54,pointsize=7.5, useDingbats = F)
par(mar = c(1, 1, 1, 1)) 

# Plot the dendrogram horizontally
plot_horiz.dendrogram(tree1, horiz = TRUE, axes = FALSE, xlim = c(0.8, 0.2), mar = c(0.01, 0.01, 0.01, 0.01))

# Add leaf-colored bar and shift it left
colored_bars(colors = my_color4(12), tree1, 
             rowLabels = ' ', 
             horiz = TRUE, 
             y_shift = -0.236)

# Draw 3 vertical colored bars manually
# Bar 1: color = color 11
rect(xleft = 0.20, xright = 0.22, 
     ybottom = 7, 
     ytop    = 12, 
     col = my_color4(12)[11], border = NA)

# Bar 2: color = color 5
rect(xleft = 0.20, xright = 0.22, 
     ybottom = 5, 
     ytop    = 6, 
     col = my_color4(12)[5], border = NA)

# Bar 3: color = color 4
rect(xleft = 0.20, xright = 0.22, 
     ybottom = 1, 
     ytop    = 4, 
     col = my_color4(12)[4], border = NA)

# Label groups
par(xpd = NA)
text(0.18, 2.5, 'Group 1', srt = 90, cex = 1)
text(0.18, 5.5, 'Group 2', srt = 90, cex = 1)
text(0.18, 9.5, 'Group 3', srt = 90, cex = 1)
par(xpd = FALSE)
dev.off()

### ================================================
### Elevation density figure 3 B
### ================================================

# Define spatial extent of the hotspot area
hotspot <- extent(361000, 880000, 2853000, 3498000)

# Mask elevation raster using county boundary
ele_focal <- mask(ele, Hds_county_outline)

# Crop to hotspot extent
ele_focal <- crop(ele_focal, hotspot)

# Resample elevation raster to match resolution and extent of beta_cata raster
ele_cal <- terra::resample(ele, beta_cata, method = 'max')

# Apply county mask to resampled elevation raster
ele_cal <- mask(ele_cal, Hds_county_outline)

# Extract elevation values for each beta diversity category
resultRaster1 <- ele_cal[beta_cata == 1]
resultRaster2 <- ele_cal[beta_cata == 2]
resultRaster3 <- ele_cal[beta_cata == 3]
resultRaster4 <- ele_cal[beta_cata == 4]
resultRaster5 <- ele_cal[beta_cata == 5]
resultRaster6 <- ele_cal[beta_cata == 6]
resultRaster7 <- ele_cal[beta_cata == 7]
resultRaster8 <- ele_cal[beta_cata == 8]
resultRaster9 <- ele_cal[beta_cata == 9]
resultRaster10 <- ele_cal[beta_cata == 10]
resultRaster11 <- ele_cal[beta_cata == 11]
resultRaster12 <- ele_cal[beta_cata == 12]

# Group beta categories into 3 major groups for plotting
resultgroup1 <- rbind(resultRaster1, resultRaster2, resultRaster3, resultRaster4)                     # Group 1
resultgroup2 <- rbind(resultRaster5, resultRaster6)                                                   # Group 2
resultgroup3 <- rbind(resultRaster7, resultRaster8, resultRaster9, resultRaster10, resultRaster11, resultRaster12)  # Group 3

# Output to PNG
#png("output/beta_cata_fig3b.png", width = 4, height = 6.9, units = "cm", res = 300, pointsize = 7.5)

# Or PDF for vector output
pdf("output/beta_cata_fig3b.pdf", width = 3.95/2.54, height = 6.9/2.54, pointsize = 7.5, useDingbats = FALSE)

# Set up plotting parameters: 1 panel, custom margins
par(mfrow = c(1, 1), mar = c(2, 2, 1, 1))

# Plot density curve for elevation values in Group 1
dens <- density(resultgroup1)
plot(dens, frame = FALSE, col = my_color4(12)[4], lwd = 1.2, 
     xlim = c(1000, 6500), ylim = c(0, 0.0015), 
     main = NA, xlab = NA, ylab = NA) 

# Overlay Group 2
dens <- density(resultgroup2)
lines(dens, frame = FALSE, col = my_color4(12)[5], lwd = 1.2, 
      xlim = c(1000, 6500), ylim = c(0, 0.0015), 
      main = NA, xlab = NA, ylab = NA, add = TRUE)

# Overlay Group 3
dens <- density(resultgroup3)
lines(dens, frame = FALSE, col = my_color4(12)[11], lwd = 1.2, 
      xlim = c(1000, 6500), ylim = c(0, 0.0015), 
      main = NA, xlab = NA, ylab = NA, add = TRUE)

# Enable plotting outside plot margins
par(xpd = NA)

# Add x-axis label ("Elevation")
text(6300, -0.00015, 'Elevation\n     (m)', cex = 0.8)

# Disable out-of-bounds plotting
par(xpd = FALSE)

# Add legend showing group-to-color mapping
legend("topright", 
       legend = c("Group 1", "Group 2", "Group 3"), 
       col = c(my_color4(12)[4], my_color4(12)[5], my_color4(12)[11]), 
       cex = 0.8, box.lty = 0, lty = 1, lwd = 2)

# Add panel label "b"
par(xpd = NA)
text(1300, 0.00155, 'b', font = 1, cex = 2)
par(xpd = FALSE)

# Add y-axis label ("Density")
par(xpd = NA)
text(1800, 0.0013, 'Density', cex = 0.8)
par(xpd = FALSE)

# Close the plotting device
dev.off()

### =========================================================================
### Plot flux vs beta Spearman figure 3 c
### =========================================================================

# read in moving window result of beta diversity
res_rel_imp <- rast('./Results/res_rel_imp_average_beta30.tif')

# Define color scale range
rn = c(-1, 1)

# Define breaks and labels for the color legend
lab = c('-1' = -1,
        "-0.5" = -0.5,
        "0" = 0,
        "0.5" = 0.5,
        "1" = 1)

# Start PNG output
#png("output/betavsflux_f3c.png", width = 7.95, height = 6.9, units = "cm", res = 300, pointsize = 7.5)

# Optional PDF output (commented out)
pdf("output/betavsflux_f3c.pdf", width = 7.95/2.54, height = 6.9/2.54,pointsize=7.5, useDingbats = F)

# Set plotting parameters
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Plot correlation raster with color scale
plot(res_rel_imp %>% terra::mask(., Hds_county_outline), 
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
     labels = "c", font = 1, cex = 2)

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

### =========================================================================
### Plot NDMS with faults and rivers figure 3 d
### =========================================================================
# Load NMDS3 raster (e.g., third axis from NMDS ordination)
NMDS3_map <- raster("./File/NMDS3_map.tif")
NMDS3_map <- rast(NMDS3_map)  # Convert to SpatRaster if needed for terra workflow

# Scale the raster values to range [-1, 1] for color normalization
scaled_raster <- ((NMDS3_map - minmax(NMDS3_map)[1]) / 
                    (minmax(NMDS3_map)[2] - minmax(NMDS3_map)[1]) - 0.5) * 2

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

# Save to PNG if needed
# png("output/NDMS_fig3D.png", width = 7.95, height = 6.9, unit = "cm", res = 300, pointsize = 7.5)

# Save to publication-quality PDF
pdf("output/NDMS_fig3d.pdf", width = 7.95 / 2.54, height = 6.9 / 2.54, pointsize = 7.5, useDingbats = FALSE)

# Set plot layout and base font size
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), ps = 8, cex = 1)

# Plot the scaled NDMS raster with no legend
plot(scaled_raster %>% terra::mask(., Hds_county_outline), 
     range = rn, 
     col = NDMS_col(17),  # 17 color levels
     plg = list(title = "NDMS3"),  # internal label, not used here
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
     labels = "d", font = 1, cex = 2)

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
                   legend.lab = "NDMS Axis 3 Values",  # legend title
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

### =========================================================================
### Figure 3 combination
### =========================================================================

# figure 3 combine
# Read in PDF pages and convert to images
# Read in the images
f3a <- image_read_pdf("output/beta_cata_fig3a.pdf", density = 300)
f3ab <- image_read_pdf("output/beta_cata_fig3ab.pdf", density = 300)
f3b <- image_read_pdf("output/beta_cata_fig3b.pdf", density = 300)
f3c <- image_read_pdf("output/betavsflux_f3c.pdf", density = 300)
f3d <- image_read_pdf("output/NDMS_fig3d.pdf", density = 300)

# Resize images to same dimensions (optional, good for alignment)
# Desired output size in cm
width_cm <- 7.95
height_cm <- 6.9

# Desired resolution (DPI)
dpi <- 300

# Convert to pixels
target_width <- round(width_cm * dpi / 2.54)
target_height <- round(height_cm * dpi / 2.54)

f3a <- image_resize(f3a, paste0(target_width, "x", target_height, "!"))
f3c <- image_resize(f3c, paste0(target_width, "x", target_height, "!"))
f3d <- image_resize(f3d, paste0(target_width, "x", target_height, "!"))
f3ab <- image_resize(f3ab, paste0(round(3.95 * dpi / 2.54), "x", target_height, "!"))
f3b <- image_resize(f3b, paste0(round(3.95 * dpi / 2.54), "x", target_height, "!"))

# Step 3: Combine rows
top_row    <- image_append(c(f3a, f3ab, f3b))
bottom_row <- image_append(c(f3c, f3d))

# Step 4: Stack the two rows vertically
fig3_combined <- image_append(c(top_row, bottom_row), stack = TRUE)

# Step 5: Save final output
image_write(fig3_combined, path = "output/Figure3_combined.png", format = "png")
image_write(fig3_combined, path = "output/Figure3_combined.pdf", format = "pdf")
