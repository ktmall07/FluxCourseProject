import pandas as pd
import numpy as np
import calc_footprint_FFP as myfootprint
from collections import Counter

###########################################################################
# Site parameters

zm = 17.575       # z-d [m]
z0 = 2.75         # roughness length [m]
h = 1250.0        # boundary layer height [m]

# Still an assumption
sigmav = 1.0

# 80% footprint contour
rs = 80.0

###########################################################################
# Read AmeriFlux file

flux_file = pd.read_csv(
    "/Users/liambogucki/Desktop/Fluxcourse-Project/AMF_US-Syv_BASE-BADM_33-5/AMF_US-Syv_BASE_HH_33-5.csv",
    encoding="cp1252",
    skiprows=2
)

# IMPORTANT: actually replace missing values
flux_file = flux_file.replace(-9999, np.nan)

###########################################################################
# Functions

def air_density(TA_C):
    """
    Air density [kg m-3]
    TA_C in deg C
    """
    return 101.3e3 / 287.04 / (TA_C + 273.15)


def calc_obukhov(rho_air, ustar, TA_C, H):
    """
    Obukhov length [m]

    TA_C : deg C
    H    : sensible heat flux [W m-2]
    """

    # avoid division by zero
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

for row in flux_file.iloc[:20].itertuples(index=False):

    rows_processed += 1

    #######################################################################
    # Use gap-filled variables where possible

    ws = row.WS_PI_F_1_1_1
    ustar = row.USTAR_PI_F_1_1_1
    wd = row.WD_PI_F_1_1_1

    #######################################################################
    # Temperature / heat flux

    TA = row.TA_PI_F_1_1_1 if hasattr(row, "TA_PI_F_1_1_1") else row.TA_1_1_1

    H = row.SH_1_1_1

    #######################################################################
    # Check required inputs

    required = [ws, ustar, wd, TA, H]

    if any(pd.isna(v) for v in required):
        rows_skipped += 1
        error_counts["Missing data"] += 1
        continue

    #######################################################################
    # FFP requires ustar > 0.1

    if ustar <= 0.1:
        rows_skipped += 1
        error_counts["ustar <= 0.1"] += 1
        continue

    #######################################################################
    # Calculate Obukhov length

    rho_air = air_density(TA)

    ol = calc_obukhov(
        rho_air=rho_air,
        ustar=ustar,
        TA_C=TA,
        H=H
    )

    if np.isnan(ol):
        rows_skipped += 1
        error_counts["H near zero"] += 1
        continue

    #######################################################################
    # FFP rejects zm/L < -15.5

    if (zm / ol) <= -15.5:
        rows_skipped += 1
        error_counts["zm/L <= -15.5"] += 1
        continue

    #######################################################################
    # Run footprint model

    try:

        FFP = myfootprint.FFP(
            zm=zm,
            z0=z0,
            h=h,
            ol=ol,
            sigmav=sigmav,
            ustar=ustar,
            wind_dir=wd,
            rs=rs,
            rslayer=1,
            fig=0          # MUCH faster
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