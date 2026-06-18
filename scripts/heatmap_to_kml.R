###############################################################################
# heatmap_to_kml.R
#
# Converts a 2D CSV heatmap grid into a shaded-contour KML for Google Earth.
#
# INPUTS:
#   1. heatmap_csv    - 2D CSV: rows = y-indices, cols = x-indices, no headers.
#   2. x_coords_csv   - 1D CSV: x-coordinates in metres, no header.
#   3. y_coords_csv   - 1D CSV: y-coordinates in metres, no header.
#   4. wind_from_deg  - Meteorological "FROM" direction (0/360 = N, 90 = E).
#                       Grid +x axis points downwind (into the wind).
#   5. origin_lat/lon - WGS-84 position of x[1], y[mid].
#   6. output_kml     - Output file path (default: "heatmap_output.kml").
#   7. n_contour_levels - Number of filled contour bands (default: 10).
#   8. palette        - RColorBrewer sequential palette (default: "YlOrRd").
#   9. layer_name     - Name shown in Google Earth layers panel.
#  10. value_label    - Quantity label used in legend and popup balloons.
#
# DEPENDENCIES: only base R + RColorBrewer (no akima, no fields, no sp/sf).
#   install.packages("RColorBrewer")
#
# HOW IT WORKS (fast path):
#   contourLines() runs on the original regular x/y metre grid — no
#   resampling or interpolation. Each contour vertex is then rotated by the
#   wind direction and projected to lon/lat with a flat-earth approximation
#   (< 1 m error over typical km-scale domains). Total cost is O(vertices),
#   not O(grid²), so even large grids finish in seconds.
#
# USAGE EXAMPLE:
#   source("heatmap_to_kml.R")
#   heatmap_to_kml(
#     heatmap_csv      = "conc.csv",
#     x_coords_csv     = "x_m.csv",
#     y_coords_csv     = "y_m.csv",
#     wind_from_deg    = 270,
#     origin_lat       = 41.8781,
#     origin_lon       = -87.6298,
#     output_kml       = "plume.kml",
#     n_contour_levels = 12,
#     palette          = "YlOrRd",
#     layer_name       = "Concentration Plume",
#     value_label      = "µg/m³"
#   )
###############################################################################

heatmap_to_kml <- function(
    heatmap_csv,
    x_coords_csv,
    y_coords_csv,
    wind_from_deg     = 0,
    origin_lat,
    origin_lon,
    output_kml        = "heatmap_output.kml",
    n_contour_levels  = 10,
    palette           = "YlOrRd",
    layer_name        = "Heatmap",
    value_label       = "value"
) {

  # ── 0. Package check ─────────────────────────────────────────────────────────
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop('Please install RColorBrewer: install.packages("RColorBrewer")')
  }

  # ── 1. Read data ──────────────────────────────────────────────────────────────
  message("Reading CSV files ...")

  Z <- as.matrix(read.csv(heatmap_csv, header = FALSE))  # [ny × nx]
  y <- unlist(read.csv(x_coords_csv,   header = FALSE))  # length nx
  x <- unlist(read.csv(y_coords_csv,   header = FALSE))  # length ny

  nx <- length(x)
  ny <- length(y)

  if (ncol(Z) != nx) stop(sprintf(
    "heatmap has %d columns but x_coords has %d values.", ncol(Z), nx))
  if (nrow(Z) != ny) stop(sprintf(
    "heatmap has %d rows but y_coords has %d values.",    nrow(Z), ny))

  message(sprintf("  Grid: %d (nx) × %d (ny)", nx, ny))
  message(sprintf("  Value range: %.4g – %.4g", min(Z, na.rm = TRUE),
                  max(Z, na.rm = TRUE)))

  # ── 2. Rotation parameters ───────────────────────────────────────────────────
  # Meteorological "FROM" bearing → rotation so local +x points downwind.
  #   wind_to  = (wind_from + 180) mod 360
  #   math_ang = 90 - wind_to        (bearing → CCW-from-east radians)
  #   rot_mat  maps local (Δx, Δy) offsets → (east_m, north_m)

  #wind_to_deg <- (wind_from_deg + 180) %% 360
  wind_to_deg <- wind_from_deg
  theta_rad   <- (90 - wind_to_deg) * pi / 180

  rot_mat <- matrix(
    c( cos(theta_rad), -sin(theta_rad),
       sin(theta_rad),  cos(theta_rad)),
    nrow = 2, byrow = TRUE
  )

  # Origin: x[1] maps to origin_lon; y[mid] maps to origin_lat
 # x_offset <- x[1]
#  y_offset  <- y[ceiling(ny / 2)]

  x_offset <- x[ceiling(nx / 2)]
  y_offset <- y[1]
  # ── 3. Flat-earth scale factors ───────────────────────────────────────────────
  # 1° lat ≈ 111 320 m;  1° lon ≈ 111 320 × cos(lat) m
  meters_per_deg_lat <- 111320
  meters_per_deg_lon <- 111320 * cos(origin_lat * pi / 180)

  # Helper: rotate a batch of local (x, y) metre-offsets → (lon, lat).
  # This is called once per contour line, so it only ever touches a few
  # hundred vertices at a time — no large matrix allocations.
  vertices_to_lonlat <- function(vx, vy) {
    dx  <- vx - x_offset
    dy  <- vy - y_offset
    rot <- rot_mat %*% rbind(dx, dy)    # [2 × n], fully vectorised
    list(
      lon = origin_lon + rot[1, ] / meters_per_deg_lon,
      lat = origin_lat + rot[2, ] / meters_per_deg_lat
    )
  }

  # ── 4. Contour lines in original grid space ───────────────────────────────────
  # contourLines() works directly on the regular metre grid — no resampling.
  # It expects z[length(x), length(y)] with x varying along rows, so transpose
  # Z (which is stored [ny × nx]).

  message("Computing contour polygons ...")

  z_min  <- min(Z, na.rm = TRUE)
  z_max  <- max(Z, na.rm = TRUE)
  breaks <- seq(z_min, z_max, length.out = n_contour_levels + 1)

  clines <- contourLines(x = x, y = y, z = t(Z), levels = breaks)

  if (length(clines) == 0) {
    warning("No contour lines generated — check that the heatmap has variation.")
  }
  message(sprintf("  Contour paths found: %d", length(clines)))

  # ── 5. Colour palette ─────────────────────────────────────────────────────────
  if (palette == "auto") palette <- "YlOrRd"
  pal_cols <- RColorBrewer::brewer.pal(max(3, min(n_contour_levels, 9)), palette)
  if (length(pal_cols) < n_contour_levels) {
    pal_cols <- colorRampPalette(pal_cols)(n_contour_levels)
  } else {
    pal_cols <- pal_cols[seq_len(n_contour_levels)]
  }

  # Convert #RRGGBB → KML aabbggrr
  hex_to_kml <- function(hex_col, alpha = "bb") {
    v <- col2rgb(hex_col)
    sprintf("%s%02x%02x%02x", alpha, v["blue", 1], v["green", 1], v["red", 1])
  }

  # ── 6. Write KML ──────────────────────────────────────────────────────────────
  message(sprintf("Writing KML to '%s' ...", output_kml))

  lines <- character(0)   # collect all KML lines into a character vector,
                          # then write once — avoids repeated paste0 growth.

  # Header
  lines <- c(lines, sprintf(
    '<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>%s</name>
    <description>Wind FROM: %g° | Origin: %.6f, %.6f | %s range: %.4g – %.4g</description>
    <open>1</open>
    <Style id="noLine">
      <LineStyle><color>00ffffff</color><width>0</width></LineStyle>
    </Style>',
    layer_name, wind_from_deg, origin_lat, origin_lon,
    value_label, z_min, z_max
  ))

  # One KML Style per contour level
  for (k in seq_len(n_contour_levels)) {
    lines <- c(lines, sprintf(
      '    <Style id="lvl%d">
      <LineStyle><color>00ffffff</color><width>0</width></LineStyle>
      <PolyStyle><color>%s</color><outline>0</outline></PolyStyle>
    </Style>',
      k, hex_to_kml(pal_cols[k])
    ))
  }

  # One Placemark per contour path
  for (cl in clines) {
    k <- findInterval(cl$level, breaks, rightmost.closed = TRUE)
    k <- max(1L, min(k, n_contour_levels))

    ll   <- vertices_to_lonlat(cl$x, cl$y)
    lons <- ll$lon
    lats <- ll$lat

    # Close ring
    if (lons[1] != lons[length(lons)] || lats[1] != lats[length(lats)]) {
      lons <- c(lons, lons[1])
      lats <- c(lats, lats[1])
    }

    coord_str <- paste(sprintf("%.6f,%.6f,0", lons, lats),
                       collapse = "\n              ")

    lines <- c(lines, sprintf(
      '    <Placemark>
      <name>%s: %.4g</name>
      <styleUrl>#lvl%d</styleUrl>
      <Polygon>
        <extrude>0</extrude>
        <altitudeMode>clampToGround</altitudeMode>
        <outerBoundaryIs><LinearRing><coordinates>
              %s
        </coordinates></LinearRing></outerBoundaryIs>
      </Polygon>
    </Placemark>',
      value_label, cl$level, k, coord_str
    ))
  }

  # Legend balloon at origin
  legend_rows <- paste(sapply(seq_len(n_contour_levels), function(k) {
    sprintf(
      '<tr><td style="background:%s;width:18px;">&nbsp;</td><td>%.3g – %.3g</td></tr>',
      pal_cols[k], breaks[k], breaks[k + 1]
    )
  }), collapse = "\n")

  lines <- c(lines, sprintf(
    '    <Placemark>
      <name>Legend</name>
      <styleUrl>#noLine</styleUrl>
      <description><![CDATA[<b>%s</b><br/><table>%s</table><br/>Wind FROM: %g°]]></description>
      <Point><coordinates>%.6f,%.6f,0</coordinates></Point>
    </Placemark>',
    value_label, legend_rows, wind_from_deg, origin_lon, origin_lat
  ))

  lines <- c(lines, "  </Document>", "</kml>")

  writeLines(lines, con = output_kml)

  message(sprintf("Done. KML written to: %s", normalizePath(output_kml)))
  invisible(output_kml)
}


###############################################################################
# ── USAGE TEMPLATE ────────────────────────────────────────────────────────────
###############################################################################

# heatmap_to_kml(
#   heatmap_csv      = "conc.csv",       # 2D matrix, no header
#   x_coords_csv     = "x_m.csv",        # 1D, metres, no header
#   y_coords_csv     = "y_m.csv",        # 1D, metres, no header
#   wind_from_deg    = 270,              # Wind FROM west → plume goes east
#   origin_lat       = 41.8781,
#   origin_lon       = -87.6298,
#   output_kml       = "plume.kml",
#   n_contour_levels = 12,
#   palette          = "YlOrRd",
#   layer_name       = "Concentration Plume",
#   value_label      = "µg/m³"
# )


###############################################################################
# ── BUILT-IN SYNTHETIC DEMO ───────────────────────────────────────────────────
# Uncomment run_demo() at the bottom to verify everything works before using
# real data. Creates demo_plume.kml in the working directory.
###############################################################################

run_demo <- function() {
  message("\n=== Synthetic Gaussian plume demo ===")

  nx <- 60; ny <- 50
  x_m <- seq(0, 3000, length.out = nx)
  y_m <- seq(-1000, 1000, length.out = ny)

  Z <- outer(x_m, y_m, function(x, y) {
    sy <- 60  + 0.12 * x
    sz <- 40  + 0.08 * x
    1e6 / (2 * pi * sy * sz) * exp(-0.5 * (y / sy)^2) * exp(-x / 2000)
  })

  tmp_z   <- tempfile(fileext = ".csv")
  tmp_x   <- tempfile(fileext = ".csv")
  tmp_y   <- tempfile(fileext = ".csv")
  out_kml <- file.path(getwd(), "demo_plume.kml")

  write.table(Z,   tmp_z, sep = ",", row.names = FALSE, col.names = FALSE)
  write.table(x_m, tmp_x, sep = ",", row.names = FALSE, col.names = FALSE)
  write.table(y_m, tmp_y, sep = ",", row.names = FALSE, col.names = FALSE)

  heatmap_to_kml(
    heatmap_csv      = tmp_z,
    x_coords_csv     = tmp_x,
    y_coords_csv     = tmp_y,
    wind_from_deg    = 270,          # FROM west → plume stretches eastward
    origin_lat       = 41.8781,
    origin_lon       = -87.6298,
    output_kml       = out_kml,
    n_contour_levels = 10,
    palette          = "YlOrRd",
    layer_name       = "Gaussian Plume Demo",
    value_label      = "µg/m³"
  )

  message(sprintf("Open in Google Earth: %s", out_kml))
  invisible(out_kml)
}

# run_demo()
