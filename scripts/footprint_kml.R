source("scripts/heatmap_to_kml.R")

args <- commandArgs(trailingOnly=TRUE)

# dir containing footprint files
fp_dir <- args[1]

# footprint files
heatmap_csv <- file.path(fp_dir, "footprint_raster_fclim2d.csv")
x_coords_csv <- file.path(fp_dir, "footprint_raster_x2d.csv")
y_coords_csv <- file.path(fp_dir, "footprint_raster_y2d.csv")
output_kml <- file.path(fp_dir, "footprint.kml")

# lat, lon, and wind direction
origin_lat <- args[2] # 41.8781
origin_lon <- args[3] # -87.6298
wind_from_deg <- args[4] # 270

# set in script
n_contour_levels <- 3
palette <- "YlOrRd"

heatmap_to_kml(
  heatmap_csv = heatmap_csv, 
  x_coords_csv = x_coords_csv, 
  y_coords_csv = y_coords_csv, 
  wind_from_deg = wind_from_deg, 
  origin_lat = origin_lat, 
  origin_lon = origin_lon, 
  output_kml = output_kml, 
  n_contour_levels = n_contour_levels, 
  palette = palette,
)