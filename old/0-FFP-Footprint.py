# Importing the needed libraries
import pandas as pd
import numpy as np
import calc_footprint_FFP as myfootprint

###########################################################################
# The params needed to generate the footprints over the time period
zm = float(
    17.575
)  # Measurement height above displacement height (i.e. z-d) [m] #Calculated from site data as 18.425m
z0 = float(
    2.75
)  # Roughness length [m] - enter [None] if not known #2.75 calculated as 0.1*27.5m
umean = 0.0  # Mean wind speed at zm [ms-1] - enter [None] if not known #WS_1_1_1
h = float(
    1250.0
)  # Boundary layer height [m] #Assume 1250m boundary height between 1000 and 1500m
ol = 0.0  # Obukhov length [m] #Need to calculate in a function
sigmav = float(
    1.0
)  # Standard deviation of lateral velocity fluctuations [ms-1] #Assuming a sigmaV of 2 as between 1 and 3, google says around 2
ustar = 0.0  # Friction velocity [ms-1] #Column called USTAR_1_1_1
wind_dir = (
    0.0  # Wind direction in degrees (of 360) for rotation of the footprint #WD_1_1_1
)
rs = float(80.0)  # Percentage of source area, i.e. a value between 10% and 90%.
# Can be either a single value (e.g., "80") or an array of increasing percentage values
# (e.g., "[10:10:80]"). Expressed either in percentages ("80") or in fractions of 1 ("0.8").
# Default is [10:10:80]. Set to "NaN" for no output of percentages.

# Reading in the US Sys sites as a test of generating the footprints
#!!!!!!!!!!!!!Will need to change path depending where the file ends up
flux_file = pd.read_csv(
    "/Users/adrianacaswell/Documents/Work/CUBioMet/Fluxcourse/FluxCourseProject/data/tower/AMF_US-Syv_BASE_HH_33-5.csv",
    encoding="cp1252",
    skiprows=2,
)

# This line was pused to verify the colums that we have to pull the variabels from in the file
# for col in flux_file.columns:
# print(col)


# Functions  to calculate the air density and the subsequent obhukov length
# Need USTAR_1_1_1, air temp (kelvin), sensible heat flux, air density, pressure
# TA_1_1_1 is temp in C so need to convert to K, so + 273.15
# H_1_1_1 is sensible heat flux
# PA_1_1_1 is the pressure
# Need a function to calculate the air density


# Function to calculate the air density
# TA is temperature
def airDensity(TA):
    return 101.3e3 / 287.04 / (TA + 273.15)


# Function to calculate the Obhukov length
def calcObhukov(airDensity, USTAR, TA, H):
    return -airDensity * (USTAR) ** 3 / (0.4 * 9.81 * ((H) / (1005.0 * (TA + 273.15))))


# Adding a counter for how many rows (half hour periods) that were skipped
rows_skipped = 0
successful_footprints = 0
unique_errors = set()
for row in flux_file.iloc[50500:50501].itertuples(index=False):

    # Running the model generation
    try:
        FFP = myfootprint.FFP(
            zm=zm,
            z0=z0,
            umean=float(row.WS_1_1_1),
            h=h,
            ol=float(
                calcObhukov(
                    airDensity(row.TA_1_1_1), row.USTAR_1_1_1, row.TA_1_1_1, row.H_1_1_1
                )
            ),
            sigmav=sigmav,
            ustar=float(row.USTAR_1_1_1),
            wind_dir=float(row.WD_1_1_1) * 100,
            rs=rs,
            fig=0,
        )
        successful_footprints += 1

        timestamp = row.TIMESTAMP_START  # or TIMESTAMP_END

        for key, value in FFP.items():
            df = pd.DataFrame(value)  # convert list of dicts → DataFrame
            df.to_csv(f"{key}.csv", index=False)

    except Exception as e:
        print(e)
        rows_skipped += 1
        continue


# Printing how many rows are skipped
print("Rows skipped")
print(rows_skipped)
print("The footprints")
print(successful_footprints)
