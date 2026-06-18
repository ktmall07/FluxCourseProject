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

# set working dir
wkdir = "/Users/adrianacaswell/Documents/Work/CUBioMet/Fluxcourse/FluxCourseProject"
os.chdir(wkdir)

# path to landcover tiff
lc_path = "data/landcover/ESA_WorldCover_10m_2021_v200_N45W090_Map.tif"

# set water value associated wiht lc product
water_val = 80

# path to footprint kml
# kml_path = "data/footprint/footprints.kml"
kml_path = "/Users/adrianacaswell/Documents/Work/CUBioMet/Fluxcourse/FluxCourseProject/footprints_US-Syv_Zoe.kml"

# read in footprint kml
fp = gpd.read_file(kml_path)

# initialize list ot add results to
rows = []

# open land cover data
with rio.open(lc_path) as lc:

    nodata = lc.nodata

    for _, row in fp.iterrows():

        ts = row["Name"]
        geom = [row["geometry"]]

        # mask landcover to footprint
        out_image, out_transform = mask(lc, geom, crop=True)

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

        # append data
        rows.append({"ts": ts, "p_water": p_water})

# convert rows to df and export to df
df = pd.DataFrame(rows)
df.to_csv("data/landcover.csv", index=False)
