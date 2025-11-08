import pyreadr
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely import wkb, wkt
from libpysal.weights import DistanceBand
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

# -----------------------------
# 0) Helpers
# -----------------------------
def to_geodf_from_rds(path):
    res = pyreadr.read_r(path)
    df = next(iter(res.values()))
    # Try WKT then WKB
    try:
        geom = gpd.GeoSeries.from_wkt(df["geometry"])
    except Exception:
        geom = df["geometry"].apply(lambda x: wkb.loads(x) if pd.notna(x) else None)
        geom = gpd.GeoSeries(geom)
    gdf = gpd.GeoDataFrame(df.drop(columns=["geometry"]), geometry=geom)
    return gdf

g_corr = to_geodf_from_rds("../Data/Shapefiles//20240627_APES FINAL_metric.rds")
