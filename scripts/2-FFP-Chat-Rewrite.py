import pandas as pd
import numpy as np
import calc_footprint_FFP as myfootprint
from collections import Counter
import simplekml
import math

###########################################################################
# Site parameters

zm = 17.575       # z-d [m]
z0 = 2.75         # roughness length [m]
h = 1250.0        # boundary layer height [m]

sigmav = 2.0

rs = 90.0

# Tower location (lat/lon) — update to actual coordinates
TOWER_LAT = 46.2421   # US-Syv latitude
TOWER_LON = -89.3477  # US-Syv longitude

###########################################################################
# Helper: offset lat/lon by dx, dy in metres

def offset_latlon(lat, lon, dx_m, dy_m):
    """
    Shift a lat/lon point by dx_m (east) and dy_m (north) metres.
    Returns (new_lat, new_lon).
    """
    R = 6378137.0  # WGS-84 Earth radius [m]
    new_lat = lat + math.degrees(dy_m / R)
    new_lon = lon + math.degrees(dx_m / (R * math.cos(math.radians(lat))))
    return new_lat, new_lon

###########################################################################
# Read AmeriFlux file

flux_file = pd.read_csv(
    "/Users/liambogucki/Desktop/Fluxcourse-Project/AMF_US-Syv_BASE-BADM_33-5/AMF_US-Syv_BASE_HH_33-5.csv",
    encoding="cp1252",
    skiprows=2
)

flux_file = flux_file.replace(-9999, np.nan)

###########################################################################
# Functions

def air_density(TA_C):
    return 101.3e3 / 287.04 / (TA_C + 273.15)


def calc_obukhov(rho_air, ustar, TA_C, H):
    if abs(H) < 1.0:
        return np.nan
    TA_K = TA_C + 273.15
    return (
        -rho_air * ustar**3
        / (0.4 * 9.81 * (H / (1005.0 * TA_K)))
    )

###########################################################################
# Counters

rows_processed = 0
rows_skipped = 0
successful_footprints = 0

error_counts = Counter()
footprints = []

###########################################################################
# Loop

for row in flux_file.iloc[0:1000].itertuples(index=False):

    rows_processed += 1

    ws    = row.WS_PI_F_1_1_1
    ustar = row.USTAR_PI_F_1_1_1
    wd    = row.WD_PI_F_1_1_1
    TA    = row.TA_PI_F_1_1_1 if hasattr(row, "TA_PI_F_1_1_1") else row.TA_1_1_1
    H     = row.SH_1_1_1

    required = [ws, ustar, wd, TA, H]

    if any(pd.isna(v) for v in required):
        rows_skipped += 1
        error_counts["Missing data"] += 1
        continue

    if ustar <= 0.1:
        rows_skipped += 1
        error_counts["ustar <= 0.1"] += 1
        continue

    rho_air = air_density(TA)
    ol = calc_obukhov(rho_air=rho_air, ustar=ustar, TA_C=TA, H=H)

    if np.isnan(ol):
        rows_skipped += 1
        error_counts["H near zero"] += 1
        continue

    if (zm / ol) <= -15.5:
        rows_skipped += 1
        error_counts["zm/L <= -15.5"] += 1
        continue

    try:
        FFP = myfootprint.FFP(
            zm=zm, z0=z0, h=h, ol=ol,
            sigmav=sigmav, ustar=ustar,
            wind_dir=wd, rs=rs,
            rslayer=1, fig=0
        )

        successful_footprints += 1
        footprints.append({
            "timestamp": row.TIMESTAMP_START,
            "footprint": FFP
        })

    except Exception as e:
        rows_skipped += 1
        error_counts[str(e)] += 1

###########################################################################
# Summary

print("\n========================")
print("PROCESSING SUMMARY")
print("========================")
print(f"Rows processed:      {rows_processed}")
print(f"Successful:          {successful_footprints}")
print(f"Skipped:             {rows_skipped}")
print("\nError summary:")
for err, count in error_counts.most_common():
    print(f"{count:6d} : {err}")
print("\nSuccess rate:")
print(f"{100 * successful_footprints / rows_processed:.2f}%")

###########################################################################
# Export footprints to KML

kml = simplekml.Kml(name="Flux Footprints - US-Syv")

for fp in footprints:
    ts  = fp["timestamp"]
    FFP = fp["footprint"]

    # FFP returns contour coordinates as lists of x and y arrays.
    # Each element of FFP['xr'] / FFP['yr'] is one contour ring
    # (there may be multiple if rs is a list; here rs is scalar so
    # we get a single-element list).

    xr_list = FFP.get("xr", [])
    yr_list = FFP.get("yr", [])

    if not xr_list or xr_list[0] is None:
        continue  # no contour returned

    for xr, yr in zip(xr_list, yr_list):
        if xr is None or yr is None:
            continue

        # Convert relative (x, y) offsets [m] to absolute lat/lon
        coords = []
        for dx, dy in zip(xr, yr):
            lat, lon = offset_latlon(TOWER_LAT, TOWER_LON, dx, dy)
            coords.append((lon, lat))  # KML uses (lon, lat) order

        # Close the ring if not already closed
        if coords[0] != coords[-1]:
            coords.append(coords[0])

        pol = kml.newpolygon(name=str(ts))
        pol.outerboundaryis = coords
        pol.style.linestyle.color  = simplekml.Color.yellow
        pol.style.linestyle.width  = 2
        pol.style.polystyle.color  = simplekml.Color.changealpha("50", simplekml.Color.yellow)
        pol.description = (
            f"Timestamp: {ts}\n"
            f"Contour: {rs}%"
        )

# Also mark the tower itself
tower_pt = kml.newpoint(name="Tower (US-Syv)")
tower_pt.coords = [(TOWER_LON, TOWER_LAT)]
tower_pt.style.iconstyle.color = simplekml.Color.red

kml_path = "footprints_US-Syv_0-1000.kml"
kml.save(kml_path)
print(f"\nKML saved to: {kml_path}")