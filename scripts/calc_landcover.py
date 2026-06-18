#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 09:19:04 2026

@author: adrianacaswell
"""

import os
import rasterio as rio
from rasterio.mask import mask
import numpy as np
import geopandas as gpd
import pandas as pd
import subprocess

# set working dir
wkdir = "/Users/adrianacaswell/Documents/Work/CUBioMet/Fluxcourse/FluxCourseProject"
os.chdir(wkdir)

# READ IN DATA
# tower lat + lon
lat = 41.8781  # 46.2421
lon = -87.6298  # -89.3459

# path to landcover tiff
lc_path = "data/landcover/ESA_WorldCover_10m_2021_v200_N45W090_Map.tif"

# set water value associated wiht lc product
water_val = 80

# path to footprint dir
fps_dir = "data/footprint"

# get dirs for each footprint timeframe
fp_dirs = [
    os.path.join(fps_dir, d)
    for d in os.listdir(fps_dir)
    if os.path.isdir(os.path.join(fps_dir, d))
]


# initialize list ot add results to
rows = []

fp_dir = fp_dirs[0]

for fp_dir in fp_dirs:

    ts = ""  # time step associated with Ameriflux data
    wdir = 270  # wind direction for the timestep

    # GET GDF FOR FOOTPRINT
    # command to run footprint_kml.R which is a wrapper for heatmap_to_R.kml
    cmd = [
        "/usr/local/bin/RScript",
        "scripts/footprint_kml.R",
        fp_dir,
        str(lat),
        str(lon),
        str(wdir),
    ]

    # create footprint.kml in fp_dir
    subprocess.run(cmd, capture_output=True, text=True)

    # read in kml file
    path_kml = f"{fp_dir}/footprint.kml"
    fp = gpd.read_file(path_kml, driver="kml")

    # MASK LANDCOVER BY FOOTPRINT
    # get footprint geometries
    geom = fp.geometry.values

    # open landcover raster
    with rio.open(lc_path) as lc:

        nodata = lc.nodata  # no data value associated with lc

        # mask
        out_image, out_transform = mask(lc, geom, crop=True)

    # CALCULATE PERCENT WATER
    # get array of raster
    arr = out_image[0]

    # mask out invalid cells
    if nodata is not None:
        arr_masked = np.ma.masked_equal(arr, nodata)

    # count number of valid cells (total area)
    n_total = arr_masked.count()

    # count number of water cells (water area)
    n_water = np.sum(arr_masked == water_val)

    # calculate percent water
    p_water = n_water / n_total

    # COLLECT DATA
    rows.append({"ts": ts, "p_water": p_water})


# convert rows to df and export to df
df = pd.DataFrame(rows)
df.to_csv("data/water.csv", index=False)


# QUICK VIS
# sample data:
# path_kml = "/Users/adrianacaswell/Documents/Work/CUBioMet/Fluxcourse/FluxCourseProject/footprints/kml/plume_30.kml"
# path_lc = "/Users/adrianacaswell/Downloads/WORLDCOVER 2/ESA_WORLDCOVER_10M_2021_V200/MAP/test_lc/ESA_WorldCover_10m_2021_v200_N45W093_Map.tif"
# dir_fp = "/Users/adrianacaswell/Documents/Work/CUBioMet/Fluxcourse/Footprint-Desai/Example 2D footprint"


# from rasterio.plot import show
# import matplotlib.pyplot as plt
# lc = rio.open(path_lc)
# fig, ax = plt.subplots(figsize = (5, 5))
# show(lc)
# fp.plot(ax=ax, facecolor = 'red', edgecolor = 'red')
